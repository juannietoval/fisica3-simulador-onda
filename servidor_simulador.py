"""
UNIVERSIDAD TECNOLÓGICA DE PEREIRA - DEPARTAMENTO DE MATEMÁTICAS
FÍSICA III - PROYECTO 10: SIMULACIÓN COMPUTACIONAL DE LA ECUACIÓN DE ONDA
========================================================================
Servidor de cálculo numérico FDTD en memoria con transmisión directa por JSON.
No genera archivos intermedios en disco duro.
"""

import json
import webbrowser
from http.server import HTTPServer, SimpleHTTPRequestHandler
import numpy as np
import os
import sys

def resolver_fdtd_en_memoria(params):
    """
    Ejecuta el método FDTD centrado en segundo orden O(dx^2, dt^2)
    estrictamente en memoria RAM y retorna los resultados en un diccionario.
    """
    # 1. Parámetros físicos
    L = float(params.get("L", 1.0))
    T = float(params.get("T", 40.0))
    mu = float(params.get("mu", 0.002))
    v = float(np.sqrt(T / mu))

    # 2. Discretización espacial
    Nx = int(params.get("Nx", 201))
    x = np.linspace(0, L, Nx)
    dx = float(x[1] - x[0])

    # 3. Discretización temporal y número de Courant (CFL)
    C = float(params.get("C", 1.0))
    dt = float(C * dx / v)
    t_max = float(params.get("t_max", 0.018))
    Nt = int(t_max / dt) + 1
    t = np.linspace(0, t_max, Nt)

    # 4. Condición inicial (t = 0)
    x0 = float(params.get("x0", 0.3 * L))
    sigma = float(params.get("sigma", 0.05 * L))
    A = float(params.get("A", 0.02))
    tipo_pulso = params.get("pulso", "gaussiano")

    y0 = np.zeros(Nx)
    if tipo_pulso == "gaussiano":
        y0 = A * np.exp(-((x - x0) ** 2) / (2.0 * (sigma ** 2)))
    elif tipo_pulso == "triangular":
        mask_izq = x <= x0
        mask_der = x > x0
        y0[mask_izq] = A * (x[mask_izq] / x0)
        y0[mask_der] = A * ((L - x[mask_der]) / (L - x0))
    elif tipo_pulso == "dos_pulsos":
        # Superposición de dos pulsos colisionantes
        y0 = (A * np.exp(-((x - 0.25 * L) ** 2) / (2.0 * (sigma ** 2))) +
              A * np.exp(-((x - 0.75 * L) ** 2) / (2.0 * (sigma ** 2))))

    bc_izq = params.get("bc_izq", "fijo")
    bc_der = params.get("bc_der", "fijo")

    def aplicar_fronteras(vec):
        if bc_izq == "fijo":
            vec[0] = 0.0
        else:
            vec[0] = vec[1]  # Extremo libre (Neumann)

        if bc_der == "fijo":
            vec[-1] = 0.0
        else:
            vec[-1] = vec[-2]  # Extremo libre (Neumann)
        return vec

    # Paso 0: t = 0
    Y = np.zeros((Nt, Nx))
    Y[0, :] = aplicar_fronteras(y0.copy())

    # Paso 1: n = 1 (arranque analítico para velocidad inicial cero)
    C2 = C ** 2
    Y[1, 1:-1] = (1.0 - C2) * Y[0, 1:-1] + 0.5 * C2 * (Y[0, 2:] + Y[0, :-2])
    Y[1, :] = aplicar_fronteras(Y[1, :])

    # Bucle temporal iterativo FDTD de segundo orden
    for n in range(1, Nt - 1):
        Y[n + 1, 1:-1] = (
            2.0 * (1.0 - C2) * Y[n, 1:-1]
            - Y[n - 1, 1:-1]
            + C2 * (Y[n, 2:] + Y[n, :-2])
        )
        Y[n + 1, :] = aplicar_fronteras(Y[n + 1, :])

    # 5. Validación analítica con d'Alembert
    idx_val = min(int(Nt * 0.45), Nt - 1)
    t_val = float(t[idx_val])
    y_analitico = np.zeros(Nx)

    if bc_izq == "fijo" and bc_der == "fijo" and tipo_pulso == "gaussiano":
        for k in range(-4, 5):
            desp = 2.0 * k * L
            f_der = A * np.exp(-((x - v * t_val - desp - x0) ** 2) / (2.0 * (sigma ** 2)))
            f_izq = A * np.exp(-((-x - v * t_val - desp - x0) ** 2) / (2.0 * (sigma ** 2)))
            f_der_reg = A * np.exp(-((x + v * t_val - desp - x0) ** 2) / (2.0 * (sigma ** 2)))
            f_izq_reg = A * np.exp(-((-x + v * t_val - desp - x0) ** 2) / (2.0 * (sigma ** 2)))
            y_analitico += 0.5 * (f_der - f_izq + f_der_reg - f_izq_reg)
        y_analitico[0] = 0.0
        y_analitico[-1] = 0.0

    diferencia = Y[idx_val, :] - y_analitico
    rmse = float(np.sqrt(np.mean(diferencia ** 2)))

    # 6. Muestreo de cuadros temporales para transmisión fluida
    total_cuadros = 120
    salto = max(1, Nt // total_cuadros)
    indices_frames = list(range(0, Nt, salto))
    
    cuadros = [Y[i, :].tolist() for i in indices_frames]
    tiempos = [float(t[i]) for i in indices_frames]

    return {
        "x": x.tolist(),
        "tiempos": tiempos,
        "cuadros": cuadros,
        "parametros": {
            "L": L, "T": T, "mu": mu, "v": round(v, 2),
            "C": round(C, 4), "dt": dt, "dx": dx, "Nt": Nt, "Nx": Nx
        },
        "validacion": {
            "t_eval": t_val,
            "y_num": Y[idx_val, :].tolist(),
            "y_anal": y_analitico.tolist(),
            "rmse_mm": float(rmse * 1e3)
        }
    }


class ManejadorSimulacion(SimpleHTTPRequestHandler):
    """Atiende peticiones HTTP estáticas y de la API JSON."""
    def __init__(self, *args, **kwargs):
        directorio_base = os.path.dirname(os.path.abspath(__file__))
        super().__init__(*args, directory=directorio_base, **kwargs)

    def do_POST(self):
        if self.path == "/api/simular":
            longitud = int(self.headers.get("Content-Length", 0))
            cuerpo = self.rfile.read(longitud)
            params = json.loads(cuerpo.decode("utf-8")) if longitud > 0 else {}

            resultado = resolver_fdtd_en_memoria(params)

            respuesta_bytes = json.dumps(resultado).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(respuesta_bytes)))
            self.send_header("Access-Control-Allow-Origin", "*")
            self.end_headers()
            self.wfile.write(respuesta_bytes)
        else:
            self.send_error(404, "Ruta no encontrada")

    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        super().end_headers()


def main():
    puerto = 8000
    directorio_actual = os.path.dirname(os.path.abspath(__file__))
    os.chdir(directorio_actual)

    servidor = HTTPServer(("localhost", puerto), ManejadorSimulacion)
    url = f"http://localhost:{puerto}/simulador_onda_web.html"

    print("=" * 80)
    print("UNIVERSIDAD TECNOLÓGICA DE PEREIRA - FÍSICA III (PROYECTO 10)")
    print("Servidor FDTD en Memoria Iniciado")
    print("=" * 80)
    print(f"Dirección local: http://localhost:{puerto}/")
    print(f"Abriendo simulador interactivo en navegador: {url}")
    print("Modo de transmisión: JSON en memoria RAM (sin archivos temporales en disco)")
    print("Presiona Ctrl + C para detener el servidor")
    print("=" * 80)

    webbrowser.open(url)
    try:
        servidor.serve_forever()
    except KeyboardInterrupt:
        print("\nServidor cerrado correctamente.")
        servidor.server_close()


if __name__ == "__main__":
    main()
