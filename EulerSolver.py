# Main routine to call and run Euler solver

import numpy as np
import matplotlib.pyplot as plt
import initialize
import RiemannSolver as riemann
import TimeStepping as tm

def lin_reconstruction(u, N):
    uL = np.zeros(N)
    uR = np.zeros(N)
    omega = 0.0

    del_uL = u[1:N-2] - u[0:N-3]
    del_uR = u[2:N-1] - u[1:N-2]
    del_i = 0.5*(1+omega)*del_uL + 0.5*(1-omega)*del_uR

    uL[1:N-2] = u[1:N-2] - 0.5*del_i
    uR[1:N-2] = u[1:N-2] + 0.5*del_i

    uL[0] = uL[1]
    uL[N-2] = uL[N-3]
    uL[N-1] = uL[N-2]

    uR[0] = uR[1]
    uR[N-2] = uR[N-3]
    uR[N-1] = uR[N-2]

    return uL, uR

def left_right_states(u):
    del_uL = u[1:N-2] - u[0:N-3]
    del_uR = u[2:N-1] - u[1:N-2]
    uL[1:N-2] = u[1:N-2] - 0.5*del_uL
    uR[1:N-2] = u[1:N-2] + 0.5*del_uL

    uL[0] = uL[1]
    uL[N-1] = uL[N-2]
    uR[0] = uR[1]
    uR[N-1] = uR[N-2]
# Call hydro
v_l = 0.1
v_r = 0.1
rho_l = 1.0
rho_r = 0.125
p_l = 1.0
p_r = 0.1
# Initial conditions
x, xmid, dx, dt, tfinal, gamma, rho, momentum, E, e_int, pressure, velocity, F, F2, F3, CFL = initialize.initi(v_l, v_r, rho_l, rho_r, p_l, p_r)
# Boundary conditions
t = 0.0
N = len(rho)

# Save initial conditions
save_rho = rho
save_velocity = velocity
save_pressure = pressure

while t<tfinal:
    rho_new = np.zeros(N)
    momentum_new = np.zeros(N)
    E_new = np.zeros(N)
    # Reconstruction to get left and right states
    rL, rR = lin_reconstruction(rho, N)
    vL, vR = lin_reconstruction(velocity, N)
    pL, pR = lin_reconstruction(pressure, N)
    # print(rL)
    # print(rR)
    p_star = pressure

    r = rho
    v = velocity
    p = pressure
    # Riemann solver
    UR_1, UR_2, UR_3, FR_1, FR_2, FR_3  = riemann.HLLE(r, rR, v, vR, p, pR, gamma, N)
    UL_1, UL_2, UL_3, FL_1, FL_2, FL_3  = riemann.HLLE(rL, r, vL, v, pL, p, gamma, N)

    U_1 = r # density
    U_2 = r * v # momentum
    U_3 = (p/(gamma-1)) + 0.5*r*(v**2) # Energy

    F_1 = U_2
    F_2 = U_2*U_1 + p
    F_3 = U_2*(U_3 + p)/U_1

    # LaxFriedrichs = 0.5*(physical_flux(uR) + physical_flux(uL) - dx/dt*(uR-uL))
    rho_new = 0.5*(FR_1 + FL_1) - (dt/dx)*(UR_1 - UL_1)
    momentum_new = 0.5*(FR_2 - FL_2) + (dt/dx)*(UR_2 - UL_2)
    E_new = 0.5*(FR_3 + FL_3) - (dt/dx)*(UR_3 - UL_3)
    # Convert back to primitives for Riemann solver
    # rL = Ul_1
    # rR = Ur_1
    # vR = Ur_2/Ur_1
    # vL = Ul_2/Ul_1
    # pR = (Ur_3 - 0.5*vR**2)*(gamma - 1)#*rR
    # pL = (Ul_3 - 0.5*vL**2)*(gamma - 1)#*rL

    # UR_1, UR_2, UR_3, FR_1, FR_2, FR_3  = riemann.HLLE(r, rR, v, vR, p, pR, gamma, N)
    # UL_1, UL_2, UL_3, FL_1, FL_2, FL_3  = riemann.HLLE(rL, r, vL, v, pL, p, gamma, N)

    # F_R, F2_R, F3_R, F_L, F2_L, F3_L  = HLLE(rL, rR, vL, vR, pL, pR, gamma, N)

    # Time-stepping
    # rho_new      = rho - (dt/dx)*(FR_1 - FL_1)
    # momentum_new = momentum - (dt/dx)*(FR_2 - FL_2)
    # E_new        = E - (dt/dx)*(FR_3 - FL_3)

    rho[:] = rho_new[:]
    momentum[:] = momentum_new[:]
    E[:] = E_new[:]
    rho[0:1] = rho[N-2:N-1]
    momentum[0:1] = momentum[N-2:N-1]
    E[0:1] = E[N-2:N-1]
    velocity = momentum/rho

    pressure = (E - 0.5*rho*velocity**2)*(gamma - 1)#*rho

    vmax = max(abs(velocity))
    if vmax <= 0.0:
        vmax = 1
    dt = CFL*dx/vmax
    t += dt
    # plt.plot(rho)
    # plt.plot(velocity)
    plt.plot(rho)
    plt.draw()
    plt.pause(0.0001)
    plt.clf()

plt.figure(1)
plt.title("Rho")
den = np.zeros(N)
den[:25] = 1.0
den[25:] = 0.125
plt.plot(save_rho)
plt.plot(rho)

plt.figure(2)
plt.title("Velocity")
plt.plot(save_velocity)
plt.plot(velocity)

plt.figure(3)
plt.title("Pressure")
plt.plot(save_pressure)
plt.plot(pressure)
plt.show()
