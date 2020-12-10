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

N = 500

x, xmid, dx, dt, tfinal, gamma, rho, momentum, E, e_int, pressure, velocity, F, F2, F3, CFL = initialize.initi(v_l, v_r, rho_l, rho_r, p_l, p_r, N)
t = 0.0
tfinal = 0.25

# Blast Wave Problem
# rho[:] = 1.0
# velocity[:] = 0.0
# # For N=500
# pressure[:50] = 1000
# pressure[50:450] = 0.1
# pressure[450:] = 100

# 123 Problem
# rho[:] = 1.0
# velocity[:250] = -2.0
# velocity[250:] = 2.0
# pressure[:250] = 0.4


# Save initial conditions
e = pressure/( rho*(gamma-1) )
energy = e + 0.5*velocity**2

save_rho = rho
save_velocity = velocity
save_pressure = pressure
save_energy = energy
nn = 0
# Sound speed
CFL=0.5
cs_2 = gamma*pressure[0]/rho[0]
while t<tfinal and nn < 2000:
    if nn%10 == 0:
        print("time: ", t)
    nn += 1
    r = rho
    v = velocity
    p = pressure
    e = p/( r*(gamma-1) )
    E = e + 0.5*v**2
    # Convert primitive variables to conservative
    ## TODO: Put variables into arrays and multiply with dot products
    q1   = r
    q2   = r*v
    q3   = r*E

    # Calculate fluxes
    F_rho   = q1*v
    F_mom1  = q2*v
    F_mom2  = p
    F_ene1  = q3*v#gamma*(q3*q2)/q1
    F_ene2  = p*v#0.5*(1-gamma)*(q2**3)/(q1**2)

    #Get left and right states of flux
    F_rhoR = np.zeros(N)
    F_rhoL = np.zeros(N)
    F_mom1R = np.zeros(N)
    F_mom1L = np.zeros(N)
    F_mom2R = np.zeros(N)
    F_mom2L = np.zeros(N)
    F_ene1R = np.zeros(N)
    F_ene1L = np.zeros(N)
    F_ene2R = np.zeros(N)
    F_ene2L = np.zeros(N)
    S_p = np.zeros(N)
    S_e = np.zeros(N)

    for i in range(1, N-1):
        if v[i] > 0.0:
            F_rhoR[i] = F_rho[i]
            F_rhoL[i] = F_rho[i-1]
            F_mom1R[i] = F_mom1[i]
            F_mom1L[i] = F_mom1[i-1]
            F_mom2R[i] = F_mom2[i]
            F_mom2L[i] = F_mom2[i-1]
            F_ene1R[i] = F_ene1[i]
            F_ene1L[i] = F_ene1[i-1]
            F_ene2R[i] = F_ene2[i]
            F_ene2L[i] = F_ene2[i-1]
        elif v[i] < 0.0:
            F_rhoR[i] = F_rho[i+1]
            F_rhoL[i] = F_rho[i]
            F_mom1R[i] = F_mom1[i+1]
            F_mom1L[i] = F_mom1[i]
            F_mom2R[i] = F_mom2[i+1]
            F_mom2L[i] = F_mom2[i]
            F_ene1R[i] = F_ene1[i+1]
            F_ene1L[i] = F_ene1[i]
            F_ene2R[i] = F_ene2[i+1]
            F_ene2L[i] = F_ene2[i]
        else:
            F_rhoR[i]  = 0.0
            F_rhoL[i]  = 0.0
            F_mom1R[i] = 0.0
            F_mom1L[i] = 0.0
            F_mom2R[i] = 0.0
            F_mom2L[i] = 0.0
            F_ene1R[i] = 0.0
            F_ene1L[i] = 0.0
            F_ene2R[i] = 0.0
            F_ene2L[i] = 0.0

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
    q3half = tm.upwind(q3, F_ene1R, F_ene1L, dt, dx)
    # Calculate source terms

    for i in range(1, N-2):
        S_p[i] = dt*(F_mom2[i+1] - F_mom2[i-1])/(2*dx)
        S_e[i] = dt*(F_ene2[i+1] - F_ene2[i-1])/(2*dx)

    S_p[0] = S_p[1]
    S_p[N-1] = S_p[N-2]
    S_e[0] = S_e[1]
    S_e[N-1] = S_e[N-2]
    # Add source terms
    q1 = q1half
    q2 = q2half - S_p
    q3 = q3half - S_e
    # Convert conservative variables to primitive
    rho = q1
    velocity = q2/q1
    # pressure = cs_2*q1
    pressure = (gamma - 1)*(q3 - 0.5*(q2**2)/q1)
    e = pressure/( rho*(gamma-1) )
    Energy = e + 0.5*velocity**2
    # pressure = (U3 - 0.5*rho*velocity**2)*(gamma - 1)
    # plt.plot(rho)
    # plt.plot(velocity)
    plt.plot(rho)
    plt.draw()
    plt.pause(0.0001)
    plt.clf()
plt.close()
# old style
# fig = plt.figure(figsize=(15,15))
fig, ax = plt.subplots(4, 1, figsize=(10,20), sharex=True)
ax[0].plot(x[:-1], rho, 'r') #row=0, col=0
ax[0].plot(x[:-1], save_rho, 'r', linestyle=":") #row=0, col=0
ax[0].set_ylabel('Density')
ax[1].plot(x[:-1], velocity, 'b') #row=0, col=0
ax[1].plot(x[:-1], save_velocity, 'b', linestyle=":") #row=0, col=0
ax[1].set_ylabel('Velocity')
ax[2].plot(x[:-1], pressure, 'g') #row=0, col=0
ax[2].plot(x[:-1], save_pressure, 'g', linestyle=":") #row=0, col=0
ax[2].set_ylabel('Pressure')
ax[3].plot(x[:-1], Energy, 'k') #row=0, col=0
ax[3].plot(x[:-1], save_energy, 'k', linestyle=":") #row=0, col=0
ax[3].set_ylabel('Energy')
plt.show()
