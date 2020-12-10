# Simple 1D convection using Finite Volume
import numpy as np

def initi (v_l, v_r, rho_l, rho_r, p_l, p_r, N):
    # Number of points
    Nx = N
    x = np.linspace(0, 1, Nx+1)
    dx = 1/Nx

    # Calculate midpoint values of x in each volume
    xmid = 0.5*(x[0:Nx] + x[1:Nx+1])
    # Set final time
    tfinal = 0.2

    # Set variables for Euler equations
    gamma = 1.66
    # Conserved Variables
    rho = np.zeros(Nx)
    momentum = np.zeros(Nx)
    E = np.zeros(Nx)

    # Primitive Variables
    e_int = np.zeros(Nx)
    pressure = np.zeros(Nx)
    velocity = np.zeros(Nx)
    F = np.zeros(Nx)
    F2 = np.zeros(Nx)
    F3 = np.zeros(Nx)

    pressure[:int(Nx/2)] = p_l
    pressure[int(Nx/2):] = p_r
    rho[:int(Nx/2)] = rho_l
    rho[int(Nx/2):] = rho_r
    velocity[:int(Nx/2)] = v_l
    velocity[int(Nx/2):] = v_r

    # Set timestep
    CFL = 0.7
    vmax = max(abs(velocity))
    if vmax <= 0.0:
        vmax = 1
    dt = CFL*dx/vmax

    return x, xmid, dx, dt, tfinal, gamma, rho, momentum, E, e_int, pressure, velocity, F, F2, F3, CFL
