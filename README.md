# Numerical Hydrodynamics
**The goal of this project is for me to apply various numerical methods to observe their pros and cons myself, as well as to implement standard testing conditions.**

## Current state of the code:
###1D Euler Equation Solver
- [x] Finite volume advection solver in 1D for the conservative, ideal Euler equations in EulerAdvection.py.
- [x] Piecewise linear reconstruction - currently uses Fromm's method.
- [x] Solves Riemman problem with HLLE solver.
- [x] Plot of Sod's shock tube problem shown in ShockTubeProblem.png

###2D Euler Equation Solver
- [x] 2D advection solver in 2D_Advection.py, using operator splitting method.
## TODO:
- [ ] Multidimensional solver for full set of Euler Equations.
- [ ] Implement HLLC solvers. I might consider adding a Roe solver as well for comparison.- [ ] Piecewise-parabolic reconstruction.
- [ ] Create better input and output (i/o) modules.
- [ ] Refactor code: Currently extremely messy, so I'll need to make it more functional.
- [ ] Addition of magnetohydrodynamics. 
