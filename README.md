# Genetic-airfoil

Optimization of airfoil geometry for a high-altitude, long endurance (HALE) UAV.

### Includes 

* Python-Xfoil interface to generate airfoil polar using Xfoil viscous analysis mode.
* Conversion from PARSEC [1] parametrization to airfoil coordinates using a Chebychev grid.
* Conversion from airfoil coordinates to PARSEC parameters through non-linear fitting via SLSQP optimization.
![plot](figs/fit.png)
* Optimization of Cl3/Cd2 endurance ratio via a genetic algorithms ## TO BE IMPLEMENTED


### References
[1] R.W. Derksen, Tim Rogalsky, “Bezier-PARSEC: An optimized aerofoil parameterization for design”, Advances in Engineering Software 41 (2010) 923-930.
