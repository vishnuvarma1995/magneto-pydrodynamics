# Numerical Hydrodynamics
**The goal of this project is for me to apply various numerical methods to observe their pros and cons myself, as well as to implement standard testing conditions.**

## Current state of the code:
-[x]I have implemented a simple finite volume advection solver for the conservative, ideal Euler equations in EulerAdvection.py. The steps leading to this were in solving the advection equation by itself in Advection.py and the isothermal equations in IsothermalAdvection.py.

## TODO:
-[]Implement option for multidimensional equation solving. I plan to do this with a simple operator splitting method. I might try other methods much later down the line. 
-[]Add approximate Riemann solvers to resolve shocks. In particular, implement HLLE and HLLC solvers. I might consider adding a Roe solver as well for comparison. 
-[]Addition of non-ideal effects and magnetohydrodynamics. 
