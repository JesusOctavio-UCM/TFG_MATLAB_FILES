# TFG: INVESTIGACIÓN MATEMÁTICA Y NUMÉRICA DE LA ECUACIÓN DE LANE-EMDEN
### Autor: Jesús Octavio Raboso
### Tutor: Uwe Brauer

Códigos empleados en la realización del TFG.

Es necesario descargar el paquete Chebfun https://www.chebfun.org/.

##### Archivos relativos al método de linealización sucesiva (SLM) para resolver la ecuacion de Lane-Emden:
- slm.m (implementación del algoritmo SLM)
  - slm_coefsA.m
  - slm_coefsB.m
  - slm_coefsR.m
  - slm_transform.m
  - cheb.m  (para computar la matriz de diferenciación de primer orden)
  - cheb2.m (para computar la matriz de diferenciación de segundo orden)
- test_slm.m (test del SLM)
- test_slm_2.m (test del SLM)

##### Archivos relativos al método de quasi-linealización sucesiva (QLM) para resolver la ecuacion de Lane-Emden:
- qlm.m (implementación del algoritmo QLM)
  - test_qlm.m (test para el algoritmo QLM)
  - test_qlm_0.m
  - test_qlm_2_5.m
  - test_qlm_cheb_coeffs.m
- qlm_first_derivative.m (implementación para la primera derivada)
  - test_qlm_first_derivative.m (test para la primera derivada)

##### Archivos para comparar distintos métodos para resolver la ecuación de Lane-Emden:
- compare_qlm_chebfun_ode45.m (compara el algoritmo QLM, la solución mediante el paquete Cehbfun y ode45)
- compare_slm_ode45.m (compara el algoritmo SLM y ode45)

##### Aachivos para generar las figuras:
- Figura 2.1: chebfun_lane_emden.m
- Figura 3.1: fdm_vs_spectral.m
- Figura 3.2: chebyshev_polynomials.m
- Figura 3.3: projection_cheb_pts.m
- Figura 3.4: aliasing.m
- Figura 3.5: gibbs_phenomenon.m
- Figura 3.6 (a): runge_phenomenon_eqispaced_pts.m
- Figura 3.6 (b): runge_phenomenon_cheb_pts.m
- Figura 3.7: example_3.m
- Figura 3.8: example_4.m
- Figura 4.1: test_qlm_0.m
- Figura 4.2: test_qlm_2_5.m
- Figura 4.3: test_qlm_cheb_coeffs.m
