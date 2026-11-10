"""
========================================================================================
UNIVERSIDAD TECNOLÓGICA DE PEREIRA - DEPARTAMENTO DE MATEMÁTICAS
FÍSICA III - EXAMEN PARCIAL 2 (PROYECTO INTEGRADOR: FENÓMENOS ONDULATORIOS)
========================================================================================

LANZADOR DEL VISOR WEB INTERACTIVO (PROYECTO 10)
Tema: Simulación computacional de la ecuación de onda unidimensional (FDTD)

Este script abre directamente en tu navegador web predeterminado la interfaz visual
interactiva con estética moderna tipo Gemini, permitiendo manipular los casos
obligatorios y crear configuraciones personalizadas en tiempo real.
========================================================================================
"""

import os
import sys
import webbrowser

def main():
    directorio_actual = os.path.dirname(os.path.abspath(__file__))
    archivo_html = os.path.join(directorio_actual, "simulador_onda_web.html")

    if not os.path.exists(archivo_html):
        print(f"Error: No se encontró el archivo '{archivo_html}'.")
        sys.exit(1)

    print("=" * 80)
    print("PROYECTO 10: LANZADOR DEL VISOR WEB INTERACTIVO")
    print("Universidad Tecnológica de Pereira - Física III (Grupo 10)")
    print("=" * 80)
    print(f"\nAbriendo el visor interactivo en tu navegador:\n{archivo_html}\n")
    print("Casos disponibles en la interfaz:")
    print("  1. Caso 1: Extremos fijos con inversión de fase (Dirichlet)")
    print("  2. Caso 2: Extremo libre sin inversión de fase (Neumann)")
    print("  3. Caso 3: Análisis de sensibilidad del parámetro de Courant CFL (0.80, 1.00, 1.05)")
    print("  4. Caso 4: Modo personalizado con controles deslizantes en tiempo real")
    print("  - Validación teórica en vivo con la solución analítica de d'Alembert (línea discontinua)")
    print("=" * 80)

    # Abrir en el navegador predeterminado
    webbrowser.open(f"file:///{os.path.abspath(archivo_html)}")

if __name__ == "__main__":
    main()
