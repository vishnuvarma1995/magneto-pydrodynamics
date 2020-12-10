# Main routine to call and solve the advection equation

import numpy as np
import matplotlib.pyplot as plt
import initialize
import RiemannSolver as riemann
import TimeStepping as tm

# Begin Program
# Initial conditions
v_l = 0.1
v_r = 0.4
rho_l = 1.0
rho_r = 0.1
p_l = 0.0
p_r = 0.0

x, xmid, dx, dt, tfinal, gamma, rho, momentum, E, e_int, pressure, velocity, F, F2, F3, CFL = initialize.initi(v_l, v_r, rho_l, rho_r, p_l, p_r)
t = 0.0
tfinal = 0.2
N = len(rho)
# Save initial conditions
save_rho = rho
save_velocity = velocity
save_pressure = pressure

nn = 0
while t<tfinal and nn < 200:
    nn += 1
    # Reconstruction to get left and right states of primitives
    r = rho
    v = velocity
    p = pressure
    e = p/(gamma - 1)

    # Conserved Quantities
    q1 = r

    # Calculate Fluxes
    v_half = np.zeros(N)
    F_rho = np.zeros(N)
    for i in range(1, N-2):
        if v[i] > 0.0:
            F_rho[i] = v[i]*q1[i-1]
        elif v[i] < 0.0:
            F_rho[i] = v[i]*q1[i]
        else:
            F_rho[i] = 0.0

    F_rho[0] = F_rho[1]
    F_rho[N-1] = F_rho[N-2]

    vmax = np.max(v)
    dt = 0.5*dx/np.abs(vmax)
    t+= dt
    # Time Stepping
    for i in range(1, N-2):
        q1[i] = q1[i] + (dt/dx)*(F_rho[i-1] - F_rho[i])
    q1[0] = q1[1]
    q1[N-1] = q1[N-2]
    # Convert conservative variables to primitive
    rho = q1
    # velocity = U2/U1

    plt.plot(rho)
    plt.draw()
    plt.pause(0.0001)
    plt.clf()

# plt.figure(1)
# plt.title("Rho")
# den = np.zeros(N)
# plt.plot(save_rho)
# plt.plot(rho)
