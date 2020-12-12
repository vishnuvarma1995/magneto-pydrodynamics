# Numerical Hydrodynamics
**The goal of this project is for me to apply various numerical methods to observe their pros and cons myself, as well as to implement standard testing conditions.**

## Current state of the code:
- [x] Finite volume advection solver in 1D for the conservative, ideal Euler equations in EulerAdvection.py.
- [x] 2D advection solver in 2D_Advection.py, using operator splitting method.

## TODO:
- [ ] Multidimensional solver for full set of Euler Equations.
- [ ] Approximate Riemann solvers to resolve shocks. In particular, implement HLLE and HLLC solvers. I might consider adding a Roe solver as well for comparison. 
- [ ] Piecewise linear reconstruction.
- [ ] Refactor code: Currently extremely messy, so I'll need to make it more functional.
- [ ] Addition of magnetohydrodynamics. 
