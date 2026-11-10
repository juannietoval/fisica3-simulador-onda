"""
========================================================================================
UNIVERSIDAD TECNOLÓGICA DE PEREIRA - DEPARTAMENTO DE MATEMÁTICAS
FÍSICA III - EXAMEN PARCIAL 2 (PROYECTO INTEGRADOR: FENÓMENOS ONDULATORIOS)
========================================================================================

PROYECTO 10: Simulación computacional de la ecuación de onda unidimensional
Grupo: 10

DESCRIPCIÓN GENERAL:
    Herramienta computacional basada en el método de Diferencias Finitas en el Dominio
    del Tiempo (FDTD) de segundo orden para modelar y simular la evolución espaciotemporal
    de perturbaciones mecánicas en una cuerda tensa bajo diversas condiciones de frontera.

MODELO FÍSICO Y MATEMÁTICO:
    Ecuación diferencial parcial hiperbólica de onda:
        ∂²y/∂t² = v² (∂²y/∂x²)
    con velocidad de propagación:
        v = sqrt(T / μ)

MÉTODO NUMÉRICO (FDTD):
    Aproximación por diferencias finitas centrales de orden O(Δx², Δt²):
        y[i, n+1] = 2*(1 - C²)*y[i, n] - y[i, n-1] + C²*(y[i+1, n] + y[i-1, n])
    donde C = v * (Δt / Δx) representa el número de Courant (condición CFL: C <= 1.0).

ESCENARIOS IMPLEMENTADOS (CUMPLIMIENTO DE EXIGENCIAS MÍNIMAS):
    1. Escenario 1: Cuerda con ambos extremos fijos (reflexión con inversión de fase).
    2. Escenario 2: Cuerda con extremo derecho libre (reflexión sin inversión de fase).
    3. Escenario 3: Análisis de sensibilidad del número de Courant (C = 1.0, 0.8, 1.05)
       demostrando exactitud, dispersión numérica e inestabilidad divergente.
    4. Validación cuantitativa: Comparación contra la solución analítica de d'Alembert
       calculando el error cuadrático medio (RMSE) y error porcentual.

========================================================================================
"""

import os
import numpy as np
import matplotlib.pyplot as plt


# ======================================================================================
# SECCIÓN 1: CONDICIONES INICIALES (FORMA Y VELOCIDAD DE LA CUERDA EN t = 0)
# ======================================================================================

def pulso_gaussiano(x, x0, sigma, A, extremos_fijos=True):
    """
    Calcula el desplazamiento inicial transversal de la cuerda en t = 0
    utilizando una función de campana gaussiana.

    Parámetros:
    -----------
    x : np.ndarray
        Vector con las posiciones de los nodos espaciales en la cuerda [m].
    x0 : float
        Posición del centro del pulso sobre la cuerda [m].
    sigma : float
        Desviación estándar o semi-ancho espacial de la campana [m].
    A : float
        Amplitud máxima de la perturbación [m].
    extremos_fijos : bool, opcional
        Si es True, fuerza a cero los desplazamientos en los extremos x=0 y x=L.
        Por defecto es True.

    Retorna:
    --------
    y0 : np.ndarray
        Vector 1D con el desplazamiento de cada nodo en el tiempo inicial [m].
    """
    y0 = A * np.exp(-((x - x0) ** 2) / (2.0 * (sigma ** 2)))
    if extremos_fijos:
        y0[0] = 0.0
        y0[-1] = 0.0
    return y0


def pulso_triangular(x, L, x_pico, A):
    """
    Genera un perfil inicial triangular continuo (modelo de cuerda pulsada).

    Parámetros:
    -----------
    x : np.ndarray
        Vector de posiciones a lo largo de la cuerda [m].
    L : float
        Longitud total de la cuerda [m].
    x_pico : float
        Posición del vértice o punto donde se levanta la cuerda [m].
    A : float
        Amplitud máxima en el vértice [m].

    Retorna:
    --------
    y0 : np.ndarray
        Vector 1D con el desplazamiento triangular inicial [m].
    """
    y0 = np.zeros_like(x)
    mascara_izq = x <= x_pico
    mascara_der = x > x_pico
    
    y0[mascara_izq] = A * (x[mascara_izq] / x_pico)
    y0[mascara_der] = A * ((L - x[mascara_der]) / (L - x_pico))
    y0[0] = 0.0
    y0[-1] = 0.0
    return y0


def velocidad_inicial(x, tipo="reposo"):
    """
    Define el perfil de velocidad transversal inicial de la cuerda en t = 0: g(x).

    Parámetros:
    -----------
    x : np.ndarray
        Vector de posiciones espaciales [m].
    tipo : str, opcional
        Tipo de condición de velocidad: "reposo" para velocidad cero en todos los puntos.
        Por defecto es "reposo".

    Retorna:
    --------
    g : np.ndarray
        Vector con la velocidad inicial en cada nodo [m/s].
    """
    if tipo == "reposo":
        return np.zeros_like(x)
    return np.zeros_like(x)


# ======================================================================================
# SECCIÓN 2: NÚCLEO DEL MÉTODO NUMÉRICO (SOLVER FDTD EN SEGUNDO ORDEN)
# ======================================================================================

def aplicar_condiciones_frontera(y, tipo_izq="fijo", tipo_der="fijo"):
    """
    Aplica las condiciones de frontera físicas en los bordes x = 0 y x = L.

    Parámetros:
    -----------
    y : np.ndarray
        Vector de desplazamientos en el instante actual [m].
    tipo_izq : str, opcional
        Condición en x = 0: "fijo" (Dirichlet, y=0) o "libre" (Neumann, dy/dx=0).
    tipo_der : str, opcional
        Condición en x = L: "fijo" o "libre".

    Retorna:
    --------
    y : np.ndarray
        Vector modificado con las condiciones de frontera aplicadas.
    """
    # Extremo izquierdo (x = 0)
    if tipo_izq == "fijo":
        y[0] = 0.0
    elif tipo_izq == "libre":
        y[0] = y[1]  # Aproximación de derivada nula dy/dx = 0

    # Extremo derecho (x = L)
    if tipo_der == "fijo":
        y[-1] = 0.0
    elif tipo_der == "libre":
        y[-1] = y[-2]  # Aproximación de derivada nula dy/dx = 0

    return y


def calcular_primer_paso(y0, g, C, dt, tipo_izq="fijo", tipo_der="fijo"):
    """
    Calcula el desplazamiento en el paso temporal n = 1 (t = dt).
    
    Justificación física:
        Para evitar el tiempo ficticio no físico y^(-1), se sustituye la derivada
        centrada de la velocidad inicial (y^1 - y^(-1))/(2*dt) = g en la ecuación
        de onda en diferencias finitas evaluada en n = 0.

    Fórmula de arranque:
        y_i^1 = y_i^0 + g_i*dt + 0.5*C²*(y_{i+1}^0 - 2*y_i^0 + y_{i-1}^0)

    Parámetros:
    -----------
    y0 : np.ndarray
        Desplazamiento inicial en t = 0 [m].
    g : np.ndarray
        Velocidad transversal inicial [m/s].
    C : float
        Número de Courant (C = v * dt / dx).
    dt : float
        Paso temporal de la discretización [s].
    tipo_izq, tipo_der : str
        Tipos de condición de frontera ("fijo" o "libre").

    Retorna:
    --------
    y1 : np.ndarray
        Desplazamiento en el paso temporal n = 1 [m].
    """
    y1 = np.zeros_like(y0)
    
    # Ecuación de recurrencia especializada en los nodos internos
    y1[1:-1] = (
        y0[1:-1]
        + g[1:-1] * dt
        + 0.5 * (C ** 2) * (y0[2:] - 2.0 * y0[1:-1] + y0[:-2])
    )
    
    # Aplicar condiciones de frontera al primer paso
    y1 = aplicar_condiciones_frontera(y1, tipo_izq, tipo_der)
    return y1


def resolver_onda_fdtd(x, t, v, y0, g, tipo_izq="fijo", tipo_der="fijo", C_forzado=None):
    """
    Ejecuta la integración numérica temporal completa de la ecuación de onda 1D
    mediante el algoritmo de Diferencias Finitas en el Dominio del Tiempo (FDTD).

    Parámetros:
    -----------
    x : np.ndarray
        Vector de discretización espacial [m].
    t : np.ndarray
        Vector de discretización temporal [s].
    v : float
        Rapidez de propagación teórica de la onda [m/s].
    y0 : np.ndarray
        Perfil de posición inicial en t = 0 [m].
    g : np.ndarray
        Perfil de velocidad inicial en t = 0 [m/s].
    tipo_izq : str, opcional
        Condición en el extremo izquierdo ("fijo" o "libre").
    tipo_der : str, opcional
        Condición en el extremo derecho ("fijo" o "libre").
    C_forzado : float, opcional
        Permite imponer un número de Courant específico para estudios de sensibilidad.

    Retorna:
    --------
    Y : np.ndarray
        Matriz bidimensional de dimensión (Nt, Nx) con la solución y(x, t) [m].
    C : float
        Número de Courant efectivo utilizado en la simulación.
    """
    Nx = len(x)
    Nt = len(t)
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    # Número de Courant: C = v * dt / dx
    if C_forzado is not None:
        C = C_forzado
    else:
        C = v * dt / dx

    # Matriz para almacenar la solución completa (filas: tiempo, columnas: espacio)
    Y = np.zeros((Nt, Nx))

    # Paso 0: Asignar condición inicial
    Y[0, :] = y0.copy()
    Y[0, :] = aplicar_condiciones_frontera(Y[0, :], tipo_izq, tipo_der)

    # Paso 1: Cálculo del primer paso temporal
    Y[1, :] = calcular_primer_paso(Y[0, :], g, C, dt, tipo_izq, tipo_der)

    # Bucle temporal iterativo para n = 1, 2, ..., Nt - 2 (calcula n + 1)
    C2 = C ** 2
    for n in range(1, Nt - 1):
        # Esquema FDTD centrado vectorial sobre todos los nodos internos
        Y[n + 1, 1:-1] = (
            2.0 * (1.0 - C2) * Y[n, 1:-1]
            - Y[n - 1, 1:-1]
            + C2 * (Y[n, 2:] + Y[n, :-2])
        )
        
        # Aplicar condiciones de frontera en el nuevo instante n + 1
        Y[n + 1, :] = aplicar_condiciones_frontera(Y[n + 1, :], tipo_izq, tipo_der)

    return Y, C


# ======================================================================================
# SECCIÓN 3: VALIDACIÓN TEÓRICA ANALÍTICA Y CÁLCULO DE ERRORES CUANTITATIVOS
# ======================================================================================

def solucion_analitica_dalembert_fija(x, t_val, v, x0, sigma, A, L, num_imagenes=5):
    """
    Calcula la solución analítica exacta de la ecuación de onda para una cuerda con
    extremos fijos mediante el principio de d'Alembert y el método de imágenes.

    Principio físico:
        Para satisfacer y(0, t) = 0 y y(L, t) = 0, la perturbación se extiende de forma
        impar y periódica con período 2*L. Cada reflexión en un extremo fijo se modela
        como una onda espejo invertida en fase (-A).

    Parámetros:
    -----------
    x : np.ndarray
        Puntos espaciales donde se evalúa la solución [m].
    t_val : float
        Instante de tiempo específico [s].
    v : float
        Rapidez de propagación [m/s].
    x0 : float
        Centro del pulso inicial [m].
    sigma : float
        Ancho del pulso inicial [m].
    A : float
        Amplitud inicial [m].
    L : float
        Longitud de la cuerda [m].
    num_imagenes : int, opcional
        Número de reflexiones/réplicas virtuales consideradas a cada lado.

    Retorna:
    --------
    y_analitico : np.ndarray
        Vector con los desplazamientos teóricos exactos en el instante t_val [m].
    """
    def f_base(coord):
        # Función pulso gaussiano base sin fronteras
        return A * np.exp(-((coord - x0) ** 2) / (2.0 * (sigma ** 2)))

    y_analitico = np.zeros_like(x)

    # Superposición de pulsos progresivos (+v) y regresivos (-v) con sus imágenes impares
    for k in range(-num_imagenes, num_imagenes + 1):
        desplazamiento = 2.0 * k * L
        
        # Onda viajera a la derecha (directa e invertida)
        y_analitico += 0.5 * f_base(x - v * t_val - desplazamiento)
        y_analitico -= 0.5 * f_base(-x - v * t_val - desplazamiento)
        
        # Onda viajera a la izquierda (directa e invertida)
        y_analitico += 0.5 * f_base(x + v * t_val - desplazamiento)
        y_analitico -= 0.5 * f_base(-x + v * t_val - desplazamiento)

    # Garantizar físicamente los extremos fijos
    y_analitico[0] = 0.0
    y_analitico[-1] = 0.0
    return y_analitico


def calcular_metricas_error(y_num, y_anal):
    """
    Calcula el Error Cuadrático Medio (RMSE) y el Error Relativo Máximo entre
    la simulación computacional y la solución teórica analítica.

    Parámetros:
    -----------
    y_num : np.ndarray
        Vector de desplazamiento numérico obtenido por FDTD [m].
    y_anal : np.ndarray
        Vector de desplazamiento analítico de d'Alembert [m].

    Retorna:
    --------
    rmse : float
        Raíz del error cuadrático medio [m].
    err_rel_pct : float
        Error relativo porcentual respecto a la amplitud máxima [%].
    """
    diferencia = y_num - y_anal
    rmse = np.sqrt(np.mean(diferencia ** 2))
    
    amplitud_max = np.max(np.abs(y_anal))
    if amplitud_max > 1e-12:
        err_rel_pct = (np.max(np.abs(diferencia)) / amplitud_max) * 100.0
    else:
        err_rel_pct = 0.0

    return rmse, err_rel_pct


# ======================================================================================
# SECCIÓN 4: GENERACIÓN DE GRÁFICAS CIENTÍFICAS DE ALTA RESOLUCIÓN
# ======================================================================================

def graficar_perfiles_instantaneos(x, t, Y, indices_tiempo, titulo, nombre_archivo):
    """
    Genera y guarda una gráfica con perfiles espaciales superpuestos en varios
    instantes de tiempo para visualizar la propagación y reflexiones.
    """
    plt.figure(figsize=(10, 5), dpi=300)
    colores = plt.cm.viridis(np.linspace(0.1, 0.9, len(indices_tiempo)))

    for idx, col in zip(indices_tiempo, colores):
        t_instante = t[idx] * 1e3  # Convertir a milisegundos
        plt.plot(
            x, Y[idx, :] * 1e3,
            label=f"t = {t_instante:.2f} ms",
            color=col,
            linewidth=1.8
        )

    plt.title(titulo, fontsize=12, fontweight="bold")
    plt.xlabel("Posición a lo largo de la cuerda x [m]", fontsize=11)
    plt.ylabel("Desplazamiento transversal y [mm]", fontsize=11)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(loc="upper right", framealpha=0.9)
    plt.tight_layout()
    plt.savefig(nombre_archivo)
    plt.close()
    print(f"Gráfica guardada: {nombre_archivo}")


def graficar_mapa_espaciotemporal(x, t, Y, titulo, nombre_archivo, v_lineas=None):
    """
    Genera un mapa de calor bidimensional en el plano espacio-tiempo (x vs t)
    que muestra la trayectoria de los paquetes de onda y sus rebotes.
    """
    plt.figure(figsize=(9, 6), dpi=300)
    X, T = np.meshgrid(x, t * 1e3)  # t en milisegundos

    val_max = np.max(np.abs(Y)) * 1e3
    if val_max == 0:
        val_max = 1.0

    mapa = plt.pcolormesh(
        X, T, Y * 1e3,
        shading="auto",
        cmap="seismic",
        vmin=-val_max,
        vmax=val_max
    )
    cbar = plt.colorbar(mapa)
    cbar.set_label("Desplazamiento y [mm]", fontsize=10)

    # Si se suministra la velocidad, trazar las trayectorias teóricas características
    if v_lineas is not None:
        t_total_ms = t[-1] * 1e3
        plt.plot([0.3, 0.3 + (v_lineas * t[-1])], [0, t_total_ms], "k--", alpha=0.5, label="Pendiente teórica (+v)")
        plt.legend(loc="upper left")

    plt.title(titulo, fontsize=12, fontweight="bold")
    plt.xlabel("Posición x [m]", fontsize=11)
    plt.ylabel("Tiempo t [ms]", fontsize=11)
    plt.tight_layout()
    plt.savefig(nombre_archivo)
    plt.close()
    print(f"Mapa espacio-temporal guardado: {nombre_archivo}")


def graficar_comparacion_analitica(x, t_val, y_num, y_anal, rmse, err_rel, nombre_archivo):
    """
    Genera una figura comparativa directa de dos paneles:
    Panel superior: Solución FDTD vs Solución d'Alembert.
    Panel inferior: Distribución espacial del error residual absoluto.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), dpi=300, sharex=True,
                                   gridspec_kw={"height_ratios": [2.5, 1.2]})

    # Panel superior: Comparación de curvas
    ax1.plot(x, y_num * 1e3, "b-", linewidth=2.0, label="FDTD Numérico (Segundo orden)")
    ax1.plot(x, y_anal * 1e3, "r--", linewidth=1.5, label="Solución Analítica de d'Alembert")
    ax1.set_title(f"Validación Numérica vs Teórica (t = {t_val*1e3:.2f} ms)", fontsize=11, fontweight="bold")
    ax1.set_ylabel("Desplazamiento y [mm]", fontsize=10)
    ax1.grid(True, linestyle="--", alpha=0.6)
    ax1.legend(loc="upper right")

    # Panel inferior: Error residual
    residual = (y_num - y_anal) * 1e3  # en mm
    ax2.plot(x, residual, "k-", linewidth=1.2, label=f"Residuo (RMSE = {rmse*1e3:.3e} mm)")
    ax2.set_xlabel("Posición x [m]", fontsize=10)
    ax2.set_ylabel("Error [mm]", fontsize=10)
    ax2.grid(True, linestyle="--", alpha=0.6)
    ax2.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(nombre_archivo)
    plt.close()
    print(f"Gráfica de validación guardada: {nombre_archivo}")


def graficar_analisis_sensibilidad(x, resultados_C, nombre_archivo):
    """
    Grafica la comparación de perfiles espaciales para diferentes valores del
    número de Courant C para evaluar estabilidad y dispersión numérica.
    """
    plt.figure(figsize=(10, 5), dpi=300)

    estilos = {"C = 1.00 (Óptimo)": ("blue", "-", 2.0),
               "C = 0.80 (Dispersivo)": ("green", "--", 1.8),
               "C = 1.05 (Inestable)": ("red", ":", 1.5)}

    for etiqueta, (y_perfil, C_val) in resultados_C.items():
        color, estilo, ancho = estilos.get(etiqueta, ("black", "-", 1.5))
        # Si el valor divergió a NaN o Infinito, acotar para graficar
        y_graficar = np.nan_to_num(y_perfil, nan=0.0, posinf=0.05, neginf=-0.05)
        plt.plot(x, y_graficar * 1e3, color=color, linestyle=estilo, linewidth=ancho, label=etiqueta)

    plt.title("Análisis de Sensibilidad: Influencia del Número de Courant (CFL)", fontsize=12, fontweight="bold")
    plt.xlabel("Posición x [m]", fontsize=11)
    plt.ylabel("Desplazamiento y [mm]", fontsize=11)
    plt.ylim(-25, 25)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(nombre_archivo)
    plt.close()
    print(f"Gráfica de sensibilidad guardada: {nombre_archivo}")


# ======================================================================================
# SECCIÓN 5: PROGRAMA PRINCIPAL (EJECUCIÓN DE LOS TRES ESCENARIOS Y ANÁLISIS)
# ======================================================================================

def main():
    print("=" * 80)
    print("PROYECTO 10: SIMULACIÓN COMPUTACIONAL DE LA ECUACIÓN DE ONDA (FÍSICA III)")
    print("=" * 80)

    # Crear carpeta para guardar los resultados gráficos junto al script
    directorio_script = os.path.dirname(os.path.abspath(__file__))
    carpeta_salida = os.path.join(directorio_script, "resultados_graficas")
    os.makedirs(carpeta_salida, exist_ok=True)

    # 1. PARÁMETROS FÍSICOS DEL SISTEMA
    L = 1.0            # Longitud de la cuerda [m]
    T = 40.0           # Tensión mecánica aplicada [N]
    mu = 0.002         # Densidad lineal de masa [kg/m]
    v = np.sqrt(T / mu)  # Rapidez teórica de propagación: v = 141.42 m/s

    # 2. MALLA ESPACIAL
    Nx = 201           # 201 nodos discretos
    x = np.linspace(0, L, Nx)
    dx = x[1] - x[0]   # dx = 5.0e-3 m = 5 mm

    # 3. CRITERIO CFL Y MALLA TEMPORAL ÓPTIMA (C = 1.0)
    C_optimo = 1.0
    dt = C_optimo * dx / v  # dt = 3.5355e-5 s
    t_max = 0.018           # Tiempo para observar viaje de ida, rebote y regreso [s]
    Nt = int(t_max / dt) + 1
    t = np.linspace(0, t_max, Nt)

    print(f"Propiedades Físicas:\n  Longitud L = {L:.2f} m | Tensión T = {T:.2f} N | Densidad mu = {mu*1e3:.2f} g/m")
    print(f"  Rapidez de onda teórica v = {v:.2f} m/s")
    print(f"Parámetros Computacionales:\n  dx = {dx*1e3:.2f} mm | dt = {dt*1e6:.2f} us | C = {C_optimo:.2f}")
    print(f"  Malla total: {Nx} puntos espaciales x {Nt} pasos temporales ({t_max*1e3:.2f} ms)\n")

    # Condición inicial: Pulso gaussiano centrado en x0 = 0.3 m, ancho sigma = 0.05 m, A = 0.02 m
    x0_pulso = 0.3
    sigma_pulso = 0.05
    A_pulso = 0.02
    y0 = pulso_gaussiano(x, x0=x0_pulso, sigma=sigma_pulso, A=A_pulso, extremos_fijos=True)
    g0 = velocidad_inicial(x, tipo="reposo")

    # ----------------------------------------------------------------------------------
    # ESCENARIO 1: EXTREMOS FIJOS (DIRICHLET) - REFLEXIÓN CON INVERSIÓN DE FASE
    # ----------------------------------------------------------------------------------
    print(">>> Ejecutando Escenario 1: Cuerda con extremos fijos...")
    Y_fijo, C1 = resolver_onda_fdtd(x, t, v, y0, g0, tipo_izq="fijo", tipo_der="fijo")

    # Selección de instantes clave para graficar
    # t = 0 ms (inicio), 2.12 ms (viaje), 4.95 ms (choque), 8.48 ms (rebote invertido)
    tiempos_deseados = [0.0, 0.0025, 0.0055, 0.0090, 0.0135]
    indices_fijo = [int(td / dt) for td in tiempos_deseados if int(td / dt) < Nt]

    graficar_perfiles_instantaneos(
        x, t, Y_fijo, indices_fijo,
        titulo="Escenario 1: Evolución Temporal con Extremos Fijos (Inversión de Fase)",
        nombre_archivo=os.path.join(carpeta_salida, "escenario1_perfiles_fijos.png")
    )
    graficar_mapa_espaciotemporal(
        x, t, Y_fijo,
        titulo="Escenario 1: Mapa Espacio-Temporal x vs t (Extremos Fijos)",
        nombre_archivo=os.path.join(carpeta_salida, "escenario1_mapa_xt.png"),
        v_lineas=v
    )

    # ----------------------------------------------------------------------------------
    # ESCENARIO 2: EXTREMO DERECHO LIBRE (NEUMANN) - REFLEXIÓN SIN INVERSIÓN
    # ----------------------------------------------------------------------------------
    print(">>> Ejecutando Escenario 2: Cuerda con extremo derecho libre...")
    Y_libre, C2 = resolver_onda_fdtd(x, t, v, y0, g0, tipo_izq="fijo", tipo_der="libre")

    graficar_perfiles_instantaneos(
        x, t, Y_libre, indices_fijo,
        titulo="Escenario 2: Evolución Temporal con Extremo Derecho Libre (Sin Inversión)",
        nombre_archivo=os.path.join(carpeta_salida, "escenario2_perfiles_libre.png")
    )
    graficar_mapa_espaciotemporal(
        x, t, Y_libre,
        titulo="Escenario 2: Mapa Espacio-Temporal x vs t (Extremo Derecho Libre)",
        nombre_archivo=os.path.join(carpeta_salida, "escenario2_mapa_xt.png")
    )

    # ----------------------------------------------------------------------------------
    # VALIDACIÓN ANALÍTICA: COMPARACIÓN CONTRA LA SOLUCIÓN DE D'ALEMBERT
    # ----------------------------------------------------------------------------------
    print(">>> Realizando Validación Analítica frente a la solución de d'Alembert...")
    t_validacion = 0.0070  # Evaluado a 7 ms (después del primer rebote)
    idx_val = int(t_validacion / dt)

    y_num_val = Y_fijo[idx_val, :]
    y_anal_val = solucion_analitica_dalembert_fija(
        x, t[idx_val], v, x0=x0_pulso, sigma=sigma_pulso, A=A_pulso, L=L
    )

    rmse_val, err_pct_val = calcular_metricas_error(y_num_val, y_anal_val)
    print(f"  Instante de prueba t = {t[idx_val]*1e3:.2f} ms:")
    print(f"  Error Cuadrático Medio (RMSE): {rmse_val*1e3:.4e} mm")
    print(f"  Error Relativo Máximo: {err_pct_val:.4f} %")

    graficar_comparacion_analitica(
        x, t[idx_val], y_num_val, y_anal_val, rmse_val, err_pct_val,
        nombre_archivo=os.path.join(carpeta_salida, "validacion_analitica_dalembert.png")
    )

    # ----------------------------------------------------------------------------------
    # ESCENARIO 3: ANÁLISIS DE SENSIBILIDAD DEL NÚMERO DE COURANT (CRITERIO CFL)
    # ----------------------------------------------------------------------------------
    print(">>> Ejecutando Escenario 3: Análisis de sensibilidad del parámetro CFL...")
    # Compararemos en un tiempo fijo t = 5.0 ms
    t_sensibilidad = 0.0050
    resultados_sensibilidad = {}

    # Caso A: C = 1.0 (Caso óptimo sin dispersión)
    Y_c1, _ = resolver_onda_fdtd(x, t, v, y0, g0, tipo_izq="fijo", tipo_der="fijo", C_forzado=1.0)
    idx_s1 = int(t_sensibilidad / dt)
    resultados_sensibilidad["C = 1.00 (Óptimo)"] = (Y_c1[idx_s1, :], 1.0)

    # Caso B: C = 0.80 (Discretización submalla: aparece dispersión numérica)
    dt_sub = 0.80 * dx / v
    Nt_sub = int(t_max / dt_sub) + 1
    t_sub = np.linspace(0, t_max, Nt_sub)
    Y_c08, _ = resolver_onda_fdtd(x, t_sub, v, y0, g0, tipo_izq="fijo", tipo_der="fijo", C_forzado=0.80)
    idx_s08 = int(t_sensibilidad / dt_sub)
    resultados_sensibilidad["C = 0.80 (Dispersivo)"] = (Y_c08[idx_s08, :], 0.80)

    # Caso C: C = 1.05 (Inestabilidad numérica divergente: C > 1.0)
    dt_sup = 1.05 * dx / v
    Nt_sup = int(t_max / dt_sup) + 1
    t_sup = np.linspace(0, t_max, Nt_sup)
    Y_c105, _ = resolver_onda_fdtd(x, t_sup, v, y0, g0, tipo_izq="fijo", tipo_der="fijo", C_forzado=1.05)
    idx_s105 = int(t_sensibilidad / dt_sup)
    resultados_sensibilidad["C = 1.05 (Inestable)"] = (Y_c105[idx_s105, :], 1.05)

    graficar_analisis_sensibilidad(
        x, resultados_sensibilidad,
        nombre_archivo=os.path.join(carpeta_salida, "escenario3_sensibilidad_cfl.png")
    )

    print("\n" + "=" * 80)
    print("SIMULACIÓN Y ANÁLISIS COMPLETADOS CON ÉXITO.")
    print(f"Todas las gráficas de alta resolución fueron generadas en: '{carpeta_salida}/'")
    print("=" * 80)


if __name__ == "__main__":
    main()
