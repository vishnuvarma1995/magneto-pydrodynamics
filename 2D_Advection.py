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

Nx = 50
Ny = 50
N = Nx
x, xmid, dx, dt, tfinal, gamma, rho, momentum, E, e_int, pressure, velocity, F, F2, F3, CFL = initialize.initi(v_l, v_r, rho_l, rho_r, p_l, p_r, N)
t = 0.0
tfinal = 1.25
dy=dx
velocityx = np.zeros((Nx, Ny))
velocityy = np.zeros((Nx, Ny))
rho = np.zeros((Nx, Ny))
pressure = np.zeros((Nx, Ny))
# v left
velocityx[:, :] = 0.1
velocityy[:, :] = 0.0
# v right
# velocityx[:, :] = 0.3
# velocityy[:, :] = 0.4
# rho left
rho[:, :] = 0.125
rho[20:30, 20:30] = 1.0
# rho right
# pressure left
pressure[:, :] = 1.0
# pressure right
# pressure[25:, :] = 0.1
# Save initial conditions
e = pressure/( rho*(gamma-1) )
energy = e + 0.5*(velocityx**2 + velocityy**2)

save_rho = rho
save_velocityx = velocityx
save_velocityy = velocityy
save_pressure = pressure
save_energy = energy
nn = 0
# Sound speed
CFL=0.5
cs_2 = gamma*pressure[0][0]/rho[0][0]

while t<tfinal and nn < 200:
    if nn%10 == 0:
        print("time: ", t)
    nn += 1
    r = rho
    vx = velocityx
    vy = velocityy
    p = pressure
    e = p/( r*(gamma-1) )
    E = e + 0.5*(vx**2 + vy**2)
    # Convert primitive variables to conservative
    ## TODO: Put variables into arrays and multiply with dot products
    q1   = r
    q2x   = r*vx
    q2y   = r*vy
    q3   = r*E

    # Calculate fluxes
    F_rhox   = q1*vx
    F_rhoy   = q1*vy
    F_mom1x  = q2x*vx # rho*u^2
    F_mom1y  = q2y*vy # rho*v^2
    F_mom1xy = q2x*vy # rho*u*v
    F_mom1yx = q2y*vx # rho*v*x
    F_mom2  = p
    F_ene1x  = q3*vx#gamma*(q3*q2)/q1
    F_ene1y  = q3*vy#gamma*(q3*q2)/q1
    F_ene2x  = p*vx#0.5*(1-gamma)*(q2**3)/(q1**2)
    F_ene2y  = p*vy#0.5*(1-gamma)*(q2**3)/(q1**2)

    #Get left and right states of flux
    F_rhoxR = np.zeros((Nx, Ny))
    F_rhoxL = np.zeros((Nx, Ny))
    F_mom1xR = np.zeros((Nx, Ny))
    F_mom1xL = np.zeros((Nx, Ny))
    F_mom1xyR = np.zeros((Nx, Ny))
    F_mom1xyL = np.zeros((Nx, Ny))
    F_ene1xR = np.zeros((Nx, Ny))
    F_ene1xL = np.zeros((Nx, Ny))
    F_ene2xR = np.zeros((Nx, Ny))
    F_ene2xL = np.zeros((Nx, Ny))

    F_mom2R = np.zeros((Nx, Ny))
    F_mom2L = np.zeros((Nx, Ny))

    F_rhoyR = np.zeros((Nx, Ny))
    F_rhoyL = np.zeros((Nx, Ny))
    F_mom1yR = np.zeros((Nx, Ny))
    F_mom1yL = np.zeros((Nx, Ny))
    F_mom1yxR = np.zeros((Nx, Ny))
    F_mom1yxL = np.zeros((Nx, Ny))
    F_ene1yR = np.zeros((Nx, Ny))
    F_ene1yL = np.zeros((Nx, Ny))
    F_ene2yR = np.zeros((Nx, Ny))
    F_ene2yL = np.zeros((Nx, Ny))

    # Calculate flux interfaces (vx-direction)
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            if vx[i][j] > 0.0:
                F_rhoxR[i][j] = F_rhox[i][j]
                F_rhoxL[i][j] = F_rhox[i-1][j]
                F_mom1xR[i][j] = F_mom1x[i][j]
                F_mom1xL[i][j] = F_mom1x[i-1][j]
                F_mom1xyR[i][j] = F_mom1xy[i][j]
                F_mom1xyL[i][j] = F_mom1xy[i-1][j]
                F_ene1xR[i][j] = F_ene1x[i][j]
                F_ene1xL[i][j] = F_ene1x[i-1][j]
                F_ene2xR[i][j] = F_ene2x[i][j]
                F_ene2xL[i][j] = F_ene2x[i-1][j]

            elif vx[i][j] < 0.0:
                F_rhoxR[i][j] = F_rhox[i+1][j]
                F_rhoxL[i][j] = F_rhox[i][j]
                F_mom1xR[i][j] = F_mom1x[i+1][j]
                F_mom1xL[i][j] = F_mom1x[i][j]
                F_mom1xyR[i][j] = F_mom1xy[i+1][j]
                F_mom1xyL[i][j] = F_mom1xy[i][j]
                F_ene1xR[i][j] = F_ene1x[i+1][j]
                F_ene1xL[i][j] = F_ene1x[i][j]
                F_ene2xR[i][j] = F_ene2x[i+1][j]
                F_ene2xL[i][j] = F_ene2x[i][j]
            else:
                F_rhoxR[i][j] = 0.0
                F_rhoxL[i][j] = 0.0
                F_mom1xR[i][j] = 0.0
                F_mom1xL[i][j] = 0.0
                F_mom1xyR[i][j] = 0.0
                F_mom1xyL[i][j] = 0.0
                F_ene1xR[i][j] = 0.0
                F_ene1xL[i][j] = 0.0
                F_ene2xR[i][j] = 0.0
                F_ene2xL[i][j] = 0.0

    # Calculate flux interfaces (vy-direction)
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            if vy[i][j] > 0.0:
                F_rhoyR[i][j] = F_rhoy[i][j]
                F_rhoyL[i][j] = F_rhoy[i][j-1]
                F_mom1yR[i][j] = F_mom1y[i][j]
                F_mom1yL[i][j] = F_mom1y[i][j-1]
                F_mom1yxR[i][j] = F_mom1yx[i][j]
                F_mom1yxL[i][j] = F_mom1yx[i][j-1]
                F_ene1yR[i][j] = F_ene1y[i][j]
                F_ene1yL[i][j] = F_ene1y[i][j-1]
                F_ene2yR[i][j] = F_ene2y[i][j]
                F_ene2yL[i][j] = F_ene2y[i][j-1]
            elif vy[i][j] < 0.0:
                F_rhoyR[i][j] = F_rhoy[i][j+1]
                F_rhoyL[i][j] = F_rhoy[i][j]
                F_mom1yR[i][j] = F_mom1y[i][j+1]
                F_mom1yL[i][j] = F_mom1y[i][j]
                F_mom1yxR[i][j] = F_mom1yx[i][j+1]
                F_mom1yxL[i][j] = F_mom1yx[i][j]
                F_ene1yR[i][j] = F_ene1y[i][j+1]
                F_ene1yL[i][j] = F_ene1y[i][j]
                F_ene2yR[i][j] = F_ene2y[i][j+1]
                F_ene2yL[i][j] = F_ene2y[i][j]
            else:
                F_rhoyR[i][j]  = 0.0
                F_rhoyL[i][j]  = 0.0
                F_mom1yR[i][j] = 0.0
                F_mom1yL[i][j] = 0.0
                F_mom1yxR[i][j] = 0.0
                F_mom1yxL[i][j] = 0.0
                F_ene1yR[i][j] = 0.0
                F_ene1yL[i][j] = 0.0
                F_ene2yR[i][j] = 0.0
                F_ene2yL[i][j] = 0.0
    # Check CFL condition
    vxmax = np.amax(abs(vx))
    vymax = np.amax(abs(vy))
    vmax = max(vxmax, vymax)
    Cs = np.sqrt(np.abs(cs_2))
    # cmax = np.amax(abs(Cs))
    dt = CFL*dx/(vmax + Cs)

    # Update in time
    t += dt
    # First advection
    q1half = np.zeros((Nx, Ny))
    q2xhalf = np.zeros((Nx, Ny))
    q2yhalf = np.zeros((Nx, Ny))
    q3half = np.zeros((Nx, Ny))
    S_p = np.zeros((Nx, Ny))
    Sx_e = np.zeros((Nx, Ny))
    Sy_e = np.zeros((Nx, Ny))

    # X-sweep
    for i in range(1, Nx-1):
        q1half[i][:] = q1[i][:] - 0.5*(dt/dx)*(F_rhoxR[i][:] - F_rhoxL[i][:])
        q2xhalf[i][:] = q2x[i][:] - 0.5*(dt/dx)*(F_mom1xR[i][:] - F_mom1xL[i][:])
        q2yhalf[i][:] = q2y[i][:] - 0.5*(dt/dx)*(F_mom1yxR[i][:] - F_mom1yxL[i][:])
        q3half[i][:] = q3[i][:] - 0.5*(dt/dx)*(F_ene1xR[i][:] - F_ene1xL[i][:]) - 0.5*(dt/dx)*(F_ene2xR[i][:] - F_ene2xL[i][:])
    # Y-sweep
    # print(q1half)
    for j in range(1, Ny-1):
        q1half[:][j] = q1half[:][j] - 0.5*(dt/dy)*(F_rhoyR[:][j] - F_rhoyL[:][j])
        q2xhalf[:][j] = q2xhalf[:][j] - 0.5*(dt/dx)*(F_mom1xyR[:][j] - F_mom1xyL[:][j])
        q2yhalf[:][j] = q2yhalf[:][j] - 0.5*(dt/dx)*(F_mom1yR[:][j] - F_mom1yL[:][j])
        q3half[:][j] = q3half[:][j] - 0.5*(dt/dx)*(F_ene1yR[:][j] - F_ene1yL[:][j]) - 0.5*(dt/dx)*(F_ene2yR[:][j] - F_ene2yL[:][j])
    # Calculate Source terms
    # print(q1half)

    # Add Source terms
    for i in range(1, Nx-1):
        S_p[i][:] = dt*(F_mom2[i+1][:] - F_mom2[i-1][:])/(2*dx)
        Sx_e[i][:] = dt*(F_ene2x[i+1][:] - F_ene2x[i-1][:])/(2*dx)
        Sy_e[i][:] = dt*(F_ene2y[i+1][:] - F_ene2y[i-1][:])/(2*dx)

    for j in range(1, Ny-1):
        S_p[:][j] = dt*(F_mom2[:][j+1] - F_mom2[:][j-1])/(2*dy)
        Sx_e[:][j] = dt*(F_ene2x[:][j+1] - F_ene2x[:][j-1])/(2*dy)
        Sy_e[:][j] = dt*(F_ene2y[:][j+1] - F_ene2y[:][j-1])/(2*dy)

    # S_p[0][:] = S_p[1][:]
    # S_p[:][0] = S_p[:][1]
    # S_p[Nx-1][:] = S_p[Nx-2][:]
    # S_p[:][Ny-1] = S_p[:][Ny-2]
    #
    # S_e[0] = S_e[1]
    # S_e[N-1] = S_e[N-2]
    # Add source terms
    q1 = q1half
    q2x = q2xhalf - S_p
    q2y = q2yhalf - S_p
    q3 = q3half - Sx_e - Sy_e

    # Y-sweep
    for j in range(1, Ny-1):
        q1[:][j] = q1[:][j] - 0.5*(dt/dy)*(F_rhoyR[:][j] - F_rhoyL[:][j])
        q2x[:][j] = q2x[:][j] - 0.5*(dt/dx)*(F_mom1xyR[:][j] - F_mom1xyL[:][j])
        q2y[:][j] = q2y[:][j] - 0.5*(dt/dx)*(F_mom1yR[:][j] - F_mom1yL[:][j])
        q3[:][j] = q3half[:][j] - 0.5*(dt/dx)*(F_ene1yR[:][j] - F_ene1yL[:][j]) - 0.5*(dt/dx)*(F_ene2yR[:][j] - F_ene2yL[:][j])
    # print(q1)

    # X-sweep
    for i in range(1, Nx-1):
        q1[i][:] = q1[i][:] - 0.5*(dt/dx)*(F_rhoxR[i][:] - F_rhoxL[i][:])
        q2x[i][:] = q2x[i][:] - 0.5*(dt/dx)*(F_mom1xR[i][:] - F_mom1xL[i][:])
        q2y[i][:] = q2y[i][:] - 0.5*(dt/dx)*(F_mom1yxR[i][:] - F_mom1yxL[i][:])
        q3[i][:] = q3[i][:] - 0.5*(dt/dx)*(F_ene1xR[i][:] - F_ene1xL[i][:])- 0.5*(dt/dx)*(F_ene2xR[i][:] - F_ene2xL[i][:])

    # print(F_rhoxR - F_rhoxL)
    # print(F_rhoyR - F_rhoyL)
    q1[0][:] = q1[1][:]
    q1[:][0] = q1[:][1]
    q1[Nx-1][:] = q1[Nx-2][:]
    q1[:][Ny-1] = q1[:][Ny-2]
    q2x[0][:] = q2x[1][:]
    q2x[:][0] = q2x[:][1]
    q2x[Nx-1][:] = q2x[Nx-2][:]
    q2x[:][Ny-1] = q2x[:][Ny-2]
    q2y[0][:] = q2y[1][:]
    q2y[:][0] = q2y[:][1]
    q2y[Nx-1][:] = q2y[Nx-2][:]
    q2y[:][Ny-1] = q2y[:][Ny-2]
    q3[0][:] = q3[1][:]
    q3[:][0] = q3[:][1]
    q3[Nx-1][:] = q3[Nx-2][:]
    q3[:][Ny-1] = q3[:][Ny-2]

    # Convert conservative variables to primitive
    # print(q1)
    rho = q1
    velocityx = q2x/q1
    velocityy = q2y/q1
    # pressure = cs_2*q1
    pressure = (gamma - 1)*(q3 - 0.5*(q2x**2+q2y**2)/q1)
    e = pressure/( rho*(gamma-1) )
    Energy = e + 0.5*(velocityx**2 + velocityy**2)
    # pressure = (U3 - 0.5*rho*velocity**2)*(gamma - 1)
    # plt.plot(rho)
    # plt.plot(velocity)
    plt.contourf(rho)
    plt.colorbar()
    plt.draw()
    plt.pause(0.0001)
    plt.clf()
# plt.close()
# old style
# fig = plt.figure(figsize=(15,15))
# fig, ax = plt.subplots(4, 1, figsize=(10,20), sharex=True)
plt.contourf(rho)
plt.colorbar()
plt.show()
