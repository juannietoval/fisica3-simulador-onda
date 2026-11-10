# Especificación Técnica del Método FDTD para la Ecuación de Onda 1D
**Asignatura:** Física III — Proyecto Integrador (Grupo 10)  
**Institución:** Universidad Tecnológica de Pereira (UTP)  
**Tema:** Simulación computacional de ondas mecánicas en una cuerda tensa

---

## 1. ¿De qué trata el método?

El método seleccionado es **Diferencias Finitas en el Dominio del Tiempo (FDTD)** centrado en segundo orden, aplicado a la ecuación de onda clásica unidimensional:

$$\frac{\partial^2 y}{\partial t^2} = v^2 \frac{\partial^2 y}{\partial x^2}$$

donde:
* $y(x, t)$ representa el desplazamiento transversal de la cuerda en la posición $x$ e instante $t$.
* $v = \sqrt{\frac{T}{\mu}}$ es la velocidad teórica de propagación de la perturbación en el medio, determinada por la tensión $T$ y la densidad lineal de masa $\mu$.

### Principio de funcionamiento
El método transforma una ecuación diferencial continua en derivadas parciales en un sistema algebraico discreto:
1. Divide la cuerda continua de longitud $L$ en $N_x$ intervalos de longitud fija $\Delta x$.
2. Divide el intervalo de observación temporal $[0, t_{\text{max}}]$ en pasos discretos de duración $\Delta t$.
3. Aproxima las derivadas segundas espacial y temporal mediante diferencias finitas centrales de orden $\mathcal{O}(\Delta x^2, \Delta t^2)$:

$$\left.\frac{\partial^2 y}{\partial t^2}\right|_i^n \approx \frac{y_i^{n+1} - 2y_i^n + y_i^{n-1}}{\Delta t^2}$$

$$\left.\frac{\partial^2 y}{\partial x^2}\right|_i^n \approx \frac{y_{i+1}^n - 2y_i^n + y_{i-1}^n}{\Delta x^2}$$

Al sustituir estas aproximaciones en la ecuación de onda se despeja el valor en el futuro inmediato $y_i^{n+1}$ en función de los valores en el presente ($y^n$) y en el pasado ($y^{n-1}$).

---

## 2. Elementos de entrada (Parámetros y Datos)

Los datos de entrada se dividen en cuatro categorías:

### A. Propiedades físicas del sistema
* **Longitud de la cuerda ($L$):** En metros ($\text{m}$).
* **Tensión mecánica ($T$):** En Newtons ($\text{N}$).
* **Densidad lineal de masa ($\mu$):** En kilogramos por metro ($\text{kg/m}$).

### B. Parámetros de discretización computacional
* **Paso espacial ($\Delta x$):** Distancia entre nodos adyacentes de la cuerda ($\text{m}$), o número de nodos $N_x = \frac{L}{\Delta x} + 1$.
* **Paso temporal ($\Delta t$):** Intervalo de avance de la simulación ($\text{s}$).
* **Tiempo total de simulación ($t_{\text{max}}$):** Duración total del análisis ($\text{s}$), con $N_t = \frac{t_{\text{max}}}{\Delta t} + 1$.

### C. Condiciones iniciales
* **Perfil de desplazamiento inicial $y(x, 0) = f(x)$:** Forma geométrica de la perturbación en $t = 0$ (ej. pulso gaussiano $A \exp\left(-\frac{(x - x_0)^2}{2\sigma^2}\right)$ o deformación triangular).
* **Perfil de velocidad inicial $\left.\frac{\partial y}{\partial t}\right|_{t=0} = g(x)$:** Velocidad transversal inicial de cada segmento (usualmente en reposo: $g(x) = 0$).

### D. Condiciones de frontera en los extremos
* **Extremo fijo (Dirichlet):** $y(0, t) = 0$ y/o $y(L, t) = 0$.
* **Extremo libre (Neumann):** $\left.\frac{\partial y}{\partial x}\right|_{x=0 \text{ ó } L} = 0$.

---

## 3. Relación Matemática entre los Datos

Los parámetros de entrada se integran a través de dos ecuaciones fundamentales:

### A. Velocidad de propagación
$$v = \sqrt{\frac{T}{\mu}}$$

### B. Parámetro de Courant (CFL)
$$C = v \frac{\Delta t}{\Delta x}$$

> **Criterio de estabilidad de Courant-Friedrichs-Lewy (CFL):**
> Para garantizar que el método sea numéricamente estable y no diverja a valores infinitos, es estrictamente obligatorio cumplir:
> $$C \le 1.0 \implies \Delta t \le \frac{\Delta x}{v}$$
> * Si $C = 1.0$: Solución numérica exacta sin dispersión artificial.
> * Si $C < 1.0$: Método estable, pero introduce dispersión numérica progresiva.
> * Si $C > 1.0$: Inestabilidad numérica inmediata (la simulación diverge).

---

## 4. Procesamiento Paso a Paso de los Datos

El algoritmo procesa los datos en cuatro etapas cronológicas:

```
[Entrada de Parámetros Físicos y de Malla]
                 │
                 ▼
[Cálculo de v = √(T/μ) y verificación de C = v·Δt/Δx ≤ 1]
                 │
                 ▼
[Paso 0: Asignación de condición inicial y(x, 0) = f(x)]
                 │
                 ▼
[Paso 1: Aplicación de condición inicial de velocidad g(x)]
                 │
                 ▼
[Bucle Temporal Iterativo para n = 1, 2, ..., Nt - 1]
  ├─ Cálculo de nodos internos mediante esquema FDTD
  └─ Imposición de condiciones de frontera (fijo/libre)
                 │
                 ▼
[Generación de Matrices y Gráficas de Salida]
```

### Detalle de las ecuaciones del algoritmo:

1. **Paso temporal $n = 0$:**
   $$y_i^0 = f(x_i), \quad \forall i \in [0, N_x - 1]$$

2. **Paso temporal inicial $n = 1$ (arranque del método):**
   Dado que $y_i^{-1}$ no existe físicamente, se utiliza la derivada central de la velocidad $\left.\frac{\partial y}{\partial t}\right|_i^0 \approx \frac{y_i^1 - y_i^{-1}}{2\Delta t} = g(x_i)$, lo que produce la fórmula especializada de arranque:
   $$y_i^1 = y_i^0 + g(x_i)\Delta t + \frac{C^2}{2} \left(y_{i+1}^0 - 2y_i^0 + y_{i-1}^0\right)$$
   Para el caso común en reposo ($g(x) = 0$):
   $$y_i^1 = (1 - C^2)y_i^0 + \frac{C^2}{2}(y_{i+1}^0 + y_{i-1}^0)$$

3. **Bucle temporal para $n \ge 1$ (puntos internos $i = 1, \dots, N_x - 2$):**
   $$y_i^{n+1} = 2(1 - C^2)y_i^n - y_i^{n-1} + C^2\left(y_{i+1}^n + y_{i-1}^n\right)$$

4. **Actualización de bordes en cada instante:**
   * Si el extremo es fijo: $y_0^{n+1} = 0$, $y_{N_x-1}^{n+1} = 0$.
   * Si el extremo $L$ es libre: $y_{N_x-1}^{n+1} = y_{N_x-2}^{n+1}$.

---

## 5. Salidas Esperadas del Sistema

El simulador produce tres tipos de resultados:

1. **Matriz de soluciones espacio-temporales:**
   Arreglo bidimensional `Y` de tamaño $(N_t, N_x)$, donde `Y[n, i]` almacena el desplazamiento $y(x_i, t_n)$.
2. **Visualizaciones científicas con unidades:**
   * **Instantáneas temporales ($y$ vs $x$):** Perfiles de la cuerda en instantes específicos ($t = 0$, $t = t_{\text{mitad}}$, $t = t_{\text{rebote}}$) que muestran la propagación y la inversión de fase en los extremos.
   * **Diagrama espacio-temporal bidimensional ($x$ vs $t$):** Mapa de calor en el que las trayectorias de las ondas se visualizan como líneas rectas con pendiente $\pm \frac{1}{v}$.
3. **Métricas de validación cuantitativa:**
   * Comparación contra la solución analítica de d'Alembert:
     $$y_{\text{analitico}}(x, t) = \frac{1}{2}\left[f(x - vt) + f(x + vt)\right]$$
   * Cálculo del Error Cuadrático Medio (RMSE) y error relativo porcentual:
     $$\text{RMSE}(t) = \sqrt{\frac{1}{N_x}\sum_{i=0}^{N_x-1} \left(y_{\text{numérico}}(x_i, t) - y_{\text{analítico}}(x_i, t)\right)^2}$$

---

## 6. Organización Modular del Código

Para mantener el código estructurado, mantenible y presentable ante el docente, se organiza en los siguientes módulos:

| Módulo | Responsabilidad | Funciones principales sugeridas |
| :--- | :--- | :--- |
| **1. Configuración (`config.py` o bloque de parámetros)** | Declarar variables físicas y discretización | $L, T, \mu, N_x, \Delta t, t_{\text{max}}$ |
| **2. Condiciones Iniciales (`initial_conditions.py`)** | Generar las funciones de perturbación | `gaussian_pulse(x, x0, sigma, A)`, `triangular_pulse(x, A)` |
| **3. Núcleo Numérico (`solver_fdtd.py`)** | Ejecutar el bucle iterativo FDTD y aplicar bordes | `solve_wave_fdtd(params, f_init, g_init, bc_type)` |
| **4. Validador Analítico (`analytical.py`)** | Calcular la solución teórica de d'Alembert | `dalembert_solution(x, t, v, f_init)` |
| **5. Postprocesamiento y Gráficas (`visualizer.py`)** | Generar gráficos de alta resolución y mapas 2D | `plot_snapshots()`, `plot_spacetime_map()`, `plot_error_analysis()` |
