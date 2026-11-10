# Simulación Computacional de la Ecuación de Onda 1D (FDTD)

**Universidad Tecnológica de Pereira (UTP)**  
**Facultad de Ciencias Básicas — Departamento de Matemáticas**  
**Asignatura:** Física III (Examen Parcial 2 — Aprendizaje Basado en Proyectos)  
**Proyecto:** 10 — Simulación computacional de la ecuación de onda  
**Grupo:** 10  

---

## 1. Descripción del Proyecto

Este proyecto implementa una herramienta computacional orientada al modelado, análisis y visualización de la propagación de ondas mecánicas transversales en una cuerda tensa unidimensional. 

Se utiliza el método de **Diferencias Finitas en el Dominio del Tiempo (FDTD)** centrado en segundo orden $\mathcal{O}(\Delta x^2, \Delta t^2)$, integrando validación cuantitativa contra la solución analítica exacta de **d'Alembert** mediante el método de imágenes y análisis de estabilidad de **Courant-Friedrichs-Lewy (CFL)**.

### Pregunta Central
> *¿Qué tan bien reproduce un modelo computacional el comportamiento de las ondas mecánicas?*

---

## 2. Modelo Físico y Matemático

La propagación de perturbaciones transversales sin amortiguamiento en una cuerda homogénea se rige por la ecuación diferencial en derivadas parciales hiperbólica:

$$\frac{\partial^2 y}{\partial t^2} = v^2 \frac{\partial^2 y}{\partial x^2}$$

donde $v$ es la rapidez teórica de fase:

$$v = \sqrt{\frac{T}{\mu}}$$

* $T$: Tensión mecánica aplicada $[\text{N}]$.
* $\mu$: Densidad lineal de masa $[\text{kg/m}]$.
* $y(x, t)$: Desplazamiento transversal $[\text{m}]$.

### Esquema Numérico FDTD (Segundo Orden)
Aproximando las derivadas espaciales y temporales mediante diferencias centrales:

$$y_i^{n+1} = 2(1 - C^2)y_i^n - y_i^{n-1} + C^2\left(y_{i+1}^n + y_{i-1}^n\right)$$

donde $C = v \frac{\Delta t}{\Delta x}$ representa el **número de Courant (CFL)**:
* $C = 1.0$: Solución numérica exacta sin dispersión artificial.
* $C < 1.0$: Esquema estable, pero introduce dispersión numérica acumulativa.
* $C > 1.0$: Inestabilidad numérica inmediata (solución divergente).

---

## 3. Estructura y Arquitectura del Proyecto

El repositorio cuenta con dos versiones complementarias e independientes:

```
fisica III python/
│
├── Fisica_III_proyecto_10.py       # Versión 1: Script científico (Matplotlib + NumPy)
├── servidor_simulador.py           # Versión 2: Servidor en memoria RAM (API JSON)
├── simulador_onda_web.html         # Versión 2: Frontend interactivo (Canvas HTML5)
├── abrir_visor_web.py              # Acceso directo al visor web
├── especificacion_metodo_fdtd.md   # Especificación técnica del método FDTD
│
└── resultados_graficas/            # Gráficas de alta resolución para el informe
    ├── escenario1_perfiles_fijos.png
    ├── escenario1_mapa_xt.png
    ├── escenario2_perfiles_libre.png
    ├── escenario2_mapa_xt.png
    ├── escenario3_sensibilidad_cfl.png
    └── validacion_analitica_dalembert.png
```

---

## 4. Instrucciones de Uso

### Requisitos
* Python 3.8 o superior.
* Librerías científicas requeridas:
  ```bash
  pip install numpy matplotlib
  ```

### Ejecución de la Versión 1 (Gráficas Científicas para el Informe)
Genera los cálculos cuantitativos y exporta las 6 figuras de alta resolución a 300 DPI:
```bash
python Fisica_III_proyecto_10.py
```

### Ejecución de la Versión 2 (Simulador Interactivo con Backend en Memoria)
Inicia el servidor HTTP ligero (sin dependencias externas de Flask/FastAPI) y transmite la simulación directamente a la memoria del navegador:
```bash
python servidor_simulador.py
```
El visor se abrirá automáticamente en:
`http://localhost:8000/simulador_onda_web.html`

---

## 5. Validación Cuantitativa

Comparación entre la solución computacional FDTD y la solución teórica de d'Alembert evaluada tras el primer ciclo de reflexión ($t = 6.97\text{ ms}$):

| Métrica | Valor Numérico | Interpretación Física |
| :--- | :--- | :--- |
| Rapidez teórica ($v$) | $141.42\text{ m/s}$ | $T = 40.0\text{ N}$, $\mu = 0.002\text{ kg/m}$. |
| Error Cuadrático Medio (RMSE) | $1.59 \times 10^{-5}\text{ mm}$ | Concordancia exacta sin error residual macroscópico. |
| Error Relativo Máximo | $0.0003\text{ \%}$ | Precisión acotada al redondeo de precisión doble. |
| Inversión en frontera fija | $\Delta\phi = \pi$ | Desfase de $180^\circ$ en extremo rígido (Dirichlet). |
| Reflexión en frontera libre | $\Delta\phi = 0$ | Conservación de amplitud y polaridad (Neumann). |

---

## 6. Créditos
Desarrollado para el Examen Parcial 2 de Física III — Universidad Tecnológica de Pereira (UTP).
