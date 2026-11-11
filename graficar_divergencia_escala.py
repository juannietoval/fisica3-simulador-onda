"""
UNIVERSIDAD TECNOLÓGICA DE PEREIRA - FÍSICA III
Demostración del Crecimiento Exponencial de la Inestabilidad Numérica (CFL > 1)
"""

import os
import numpy as np
import matplotlib.pyplot as plt

def main():
    directorio_actual = os.path.dirname(os.path.abspath(__file__))
    carpeta_salida = os.path.join(directorio_actual, "resultados_graficas")
    os.makedirs(carpeta_salida, exist_ok=True)

    # Parámetros físicos base
    L = 1.0
    T = 40.0
    mu = 0.002
    v = np.sqrt(T / mu) # 141.42 m/s

    Nx = 151
    x = np.linspace(0, L, Nx)
    dx = x[1] - x[0]

    # C = 1.05 (Inestable: supera el límite CFL)
    C = 1.05
    dt = C * dx / v
    t_max = 0.012 # 12 ms
    Nt = int(t_max / dt) + 1
    t = np.linspace(0, t_max, Nt)

    # Condición inicial: pulso gaussiano suave A = 20 mm
    A = 0.02
    x0 = 0.4
    sigma = 0.06
    y0 = A * np.exp(-((x - x0) ** 2) / (2.0 * (sigma ** 2)))
    y0[0] = 0.0
    y0[-1] = 0.0

    # Arreglo FDTD
    Y = np.zeros((Nt, Nx))
    Y[0, :] = y0

    # Paso 1
    C2 = C ** 2
    Y[1, 1:-1] = (1.0 - C2) * Y[0, 1:-1] + 0.5 * C2 * (Y[0, 2:] + Y[0, :-2])
    Y[1, 0] = 0.0
    Y[1, -1] = 0.0

    max_amplitudes = [np.max(np.abs(Y[0, :]))]
    max_amplitudes.append(np.max(np.abs(Y[1, :])))

    # Bucle temporal
    for n in range(1, Nt - 1):
        Y[n + 1, 1:-1] = (
            2.0 * (1.0 - C2) * Y[n, 1:-1]
            - Y[n - 1, 1:-1]
            + C2 * (Y[n, 2:] + Y[n, :-2])
        )
        Y[n + 1, 0] = 0.0
        Y[n + 1, -1] = 0.0
        max_amplitudes.append(np.max(np.abs(Y[n + 1, :])))

    max_amplitudes = np.array(max_amplitudes)

    # Crear Figura Científica Comparativa de 2 Paneles
    fig = plt.figure(figsize=(11, 7), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1])

    # Subplot 1: Perfiles iniciales y primeros instantes donde nace la inestabilidad
    ax1 = fig.add_subplot(gs[0, 0])
    instantes_iniciales = [0.0, 0.002, 0.004, 0.0055]
    colores1 = ["#2563eb", "#059669", "#d97706", "#dc2626"]
    for ti, col in zip(instantes_iniciales, colores1):
        idx = int(ti / dt)
        if idx < Nt:
            ax1.plot(x, Y[idx, :] * 1e3, label=f"t = {ti*1e3:.1f} ms", color=col, linewidth=1.6)
    ax1.set_title("Fase Inicial (t ≤ 5.5 ms): El pulso se deforma", fontsize=10, fontweight="bold")
    ax1.set_xlabel("Posición x [m]", fontsize=9)
    ax1.set_ylabel("Desplazamiento [mm]", fontsize=9)
    ax1.set_ylim(-35, 35)
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="upper right", fontsize=8)

    # Subplot 2: Explosión numérica en escala grande (Metros)
    ax2 = fig.add_subplot(gs[0, 1])
    instantes_explosion = [0.0065, 0.0075, 0.0085]
    colores2 = ["#7c3aed", "#e11d48", "#991b1b"]
    for ti, col in zip(instantes_explosion, colores2):
        idx = int(ti / dt)
        if idx < Nt:
            # En metros
            ax2.plot(x, Y[idx, :], label=f"t = {ti*1e3:.1f} ms (pico: {np.max(np.abs(Y[idx, :])):.1f} m)",
                     color=col, linewidth=1.4)
    ax2.set_title("Fase de Explosión: Amplitud crece a METROS", fontsize=10, fontweight="bold")
    ax2.set_xlabel("Posición x [m]", fontsize=9)
    ax2.set_ylabel("Desplazamiento [m]", fontsize=9)
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="upper left", fontsize=8)

    # Subplot 3: Curva de Crecimiento Exponencial (Semilogarítmica)
    ax3 = fig.add_subplot(gs[1, :])
    ax3.semilogy(t * 1e3, max_amplitudes * 1e3, color="#b91c1c", linewidth=2.0, label="Amplitud Máxima |y|_max (t)")
    ax3.axhline(A * 1e3, color="black", linestyle=":", label="Amplitud física inicial (20 mm)")
    ax3.set_title("Crecimiento Exponencial del Error Numérico (Escala Logarítmica)", fontsize=11, fontweight="bold")
    ax3.set_xlabel("Tiempo transcurrido [ms]", fontsize=10)
    ax3.set_ylabel("Amplitud máxima [mm] (Log)", fontsize=10)
    ax3.grid(True, which="both", linestyle="--", alpha=0.5)
    ax3.legend(loc="upper left")

    plt.tight_layout()
    ruta_guardado = os.path.join(carpeta_salida, "crecimiento_inestabilidad_escala.png")
    plt.savefig(ruta_guardado)
    plt.close()
    print(f"Gráfica guardada exitosamente en: {ruta_guardado}")

if __name__ == "__main__":
    main()
