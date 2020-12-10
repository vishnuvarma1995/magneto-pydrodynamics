# Main routine to call and solve the advection equation

import numpy as np
import matplotlib.pyplot as plt
import initialize
import RiemannSolver as riemann
import TimeStepping as tm

# Begin Program
# Initial conditions
v_l = 0.0
v_r = 0.0
rho_l = 1.0
rho_r = 0.125
p_l = 1.0
p_r = 0.1

x, xmid, dx, dt, tfinal, gamma, rho, momentum, E, e_int, pressure, velocity, F, F2, F3, CFL = initialize.initi(v_l, v_r, rho_l, rho_r, p_l, p_r)
t = 0.0
tfinal = 1.0
N = len(rho)
# Save initial conditions
save_rho = rho
save_velocity = velocity
save_pressure = pressure
nn = 0
# Sound speed
cs_2 = gamma*pressure[0]/rho[0]
while t<tfinal and nn < 200:
    nn += 1

    r = rho
    v = velocity
    p = pressure

    # Convert primitive variables to conservative
    ## TODO: Put variables into arrays and multiply with dot products
    q1   = r
    q2   = r*v

    # Calculate fluxes
    F_rho  = r*v
    F_mom1  = r*v*v
    F_mom2  = p


    #Get left and right states of flux
    F_rhoR = np.zeros(N)
    F_rhoL = np.zeros(N)
    F_mom1R = np.zeros(N)
    F_mom1L = np.zeros(N)
    F_mom2R = np.zeros(N)
    F_mom2L = np.zeros(N)
    S_p = np.zeros(N)

    for i in range(1, N-1):
        if v[i] > 0.0:
            F_rhoR[i] = F_rho[i]
            F_rhoL[i] = F_rho[i-1]
            F_mom1R[i] = F_mom1[i]
            F_mom1L[i] = F_mom1[i-1]
            F_mom2R[i] = F_mom2[i]
            F_mom2L[i] = F_mom2[i-1]
        elif v[i] < 0.0:
            F_rhoR[i] = F_rho[i+1]
            F_rhoL[i] = F_rho[i]
            F_mom1R[i] = F_mom1[i+1]
            F_mom1L[i] = F_mom1[i]
            F_mom2R[i] = F_mom2[i+1]
            F_mom2L[i] = F_mom2[i]
        else:
            F_rhoR[i] = 0.0
            F_rhoL[i] = 0.0
            F_mom1R[i] = 0.0
            F_mom1L[i] = 0.0
            F_mom2R[i] = 0.0
            F_mom2L[i] = 0.0

    # Check CFL condition
    vmax = max(abs(v))
    Cs = np.sqrt(np.abs(cs_2))
    # cmax = max(abs(Cs))
    dt = CFL*dx/(vmax + Cs)

    # Update in time
    t += dt
    # First advection
    # for i in range(1, N-2):
    #     q1half = q1 - (dt/dx)*(F_rho[i+1] - F_rho[i])
    #     q2half = q2 - (dt/dx)*(F_mom1[i+1] - F_mom1[i])
    # q1half[0] = q1half[1]
    # q1half[N-1] = q1half[N-2]
    # q2half[0] = q2half[1]
    # q2half[N-1] = q2half[N-2]
    q1half = tm.upwind(q1, F_rhoR, F_rhoL, dt, dx)
    q2half = tm.upwind(q2, F_mom1R, F_mom1L, dt, dx)
    # Calculate source terms

    for i in range(1, N-2):
        S_p[i] = cs_2*(F_mom2[i+1] - F_mom2[i-1])/(2*dx)

    S_p[0] = S_p[1]
    S_p[N-1] = S_p[N-2]
    # Add source terms
    q1 = q1half
    q2 = q2half - S_p
    # Convert conservative variables to primitive
    rho = q1
    velocity = q2/q1
    pressure = cs_2*q1
    # pressure = (U3 - 0.5*rho*velocity**2)*(gamma - 1)
    # plt.plot(rho)
    # plt.plot(velocity)
    plt.plot(rho)
    plt.draw()
    plt.pause(0.0001)
    plt.clf()
