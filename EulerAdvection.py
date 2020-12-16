# Main routine to call and solve the advection equation

import numpy as np
import matplotlib.pyplot as plt
import initialize
import RiemannSolver as riemann
import TimeStepping as tm
import Reconstruction as reconstruct

# Begin Program
# Initial conditions
# v_l = 0.0
# v_r = 0.0
# rho_l = 1.0
# rho_r = 0.125
# p_l = 1.0
# p_r = 0.1

N = 500
Nx = N
gamma = 1.66 # For gamma law pressure, 1.66 for ideal gas
dx = 1/Nx
x = np.linspace(0, 1, Nx)

# x, xmid, dx, dt, tfinal, gamma, rho, momentum, E, e_int, pressure, velocity, F, F2, F3, CFL = initialize.initi(v_l, v_r, rho_l, rho_r, p_l, p_r, N)
t = 0.0
tfinal = 0.25
rho = np.zeros(Nx)
energy = np.zeros(Nx)
pressure = np.zeros(Nx)
velocity = np.zeros(Nx)

# Sod's shock tube problem
pressure[:int(Nx/2)] = 1.0
pressure[int(Nx/2):] = 0.1
rho[:int(Nx/2)] = 1.0
rho[int(Nx/2):] = 0.125
velocity[:int(Nx/2)] = 0.0
velocity[int(Nx/2):] = 0.0

#Blast Wave Problem
# rho[:] = 1.0
# velocity[:] = 0.0
# pressure[:int(Nx/10)] = 1000
# pressure[int(Nx/10):int(9*Nx/10)] = 0.1
# pressure[int(9*Nx/10):] = 100

# 123 Problem
# rho[:] = 1.0
# velocity[:int(Nx/2)] = -2.0
# velocity[int(Nx/2):] = 2.0
# pressure[:] = 0.4

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

    # Calculate sound speed
    Cs = np.sqrt(gamma*np.array(np.abs(p))/np.array(np.abs(r)))

    # Get left and right states of flux
    F_rhoR, F_rhoL, F_mom1R, F_mom1L, F_ene1R, F_ene1L = reconstruct.HLLE_reconstruction_1D(q1,q2,q3,v,Cs,N,dx, reconstruction='lin', linearmethod='centered')

    # Check CFL condition
    dt = CFL*dx/(np.max(np.abs(v+Cs)))
    # Initialize arrays for source terms (variables treated as source terms to advection problem)
    S_p = np.zeros(N)
    S_e = np.zeros(N)
    # Update in time
    t += dt
    q1half = tm.upwind(q1, F_rhoR, F_rhoL, dt, dx)
    q2half = tm.upwind(q2, F_mom1R, F_mom1L, dt, dx)
    q3half = tm.upwind(q3, F_ene1R, F_ene1L, dt, dx)
    # Update source terms
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
    pressure = (gamma - 1)*(q3 - 0.5*(q2**2)/q1)
    e = pressure/( rho*(gamma-1) )
    Energy = e + 0.5*velocity**2

#     plt.plot(rho)
#     plt.draw()
#     plt.pause(0.0001)
#     plt.clf()
# plt.close()

fig, ax = plt.subplots(4, 1, figsize=(10,20), sharex=True)
ax[0].plot(x[:], rho, 'r')
ax[0].plot(x[:], save_rho, 'r', linestyle=":")
ax[0].set_ylabel('Density')
ax[1].plot(x[:], velocity, 'b')
ax[1].plot(x[:], save_velocity, 'b', linestyle=":")
ax[1].set_ylabel('Velocity')
ax[2].plot(x[:], pressure, 'g')
ax[2].plot(x[:], save_pressure, 'g', linestyle=":")
ax[2].set_ylabel('Pressure')
ax[3].plot(x[:], Energy, 'k')
ax[3].plot(x[:], save_energy, 'k', linestyle=":")
ax[3].set_ylabel('Energy')
plt.show()
