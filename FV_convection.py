# Simple 1D convection using Finite Volume
import numpy as np
import matplotlib.pyplot as plt

# Number of points
Nx = 50
x = np.linspace(0, 1, Nx+1)
dx = 1/Nx

# Calculate midpoint values of x in each volume
xmid = 0.5*(x[0:Nx] + x[1:Nx+1])

# Set velocity
#u = 0.8

# Set final time
tfinal = 0.2


# Set variables for Euler equations
gamma = 1.66
rho = np.zeros(Nx)
momentum = np.zeros(Nx)
eos = np.zeros(Nx)
e = np.zeros(Nx)
Energy = np.zeros(Nx)
pressure = np.zeros(Nx)
F = np.zeros(Nx)
F2 = np.zeros(Nx)
F3 = np.zeros(Nx)
velocity = np.zeros(Nx)
#velocity[:int(Nx/2)] = 0.3
#velocity[int(Nx/2):] = -0.1
#velocity = 0.75*np.exp(-((xmid - 0.5)/0.1)**2)
#pressure[:int(Nx/2)] = 100
pressure[:int(Nx/2)] = 1
pressure[int(Nx/2):] = 0.1
rho[:int(Nx/2)] = 1
rho[int(Nx/2):] = 0.125
e = pressure/((gamma - 1)*rho)
Energy = rho*e + 0.5*rho*velocity**2


# Set timestep
CFL = 0.7
dt = CFL*dx#/np.max(np.abs(velocity))


def upwind_flux(U, velocity,dt, dx, Nx, F_rho, F_mom, F_ene, pressure, Energy):
    Ur = np.roll(U, 1)
    Ul = np.roll(U, -1)
    R =np.zeros(Nx)
    R_mom =np.zeros(Nx)
    R_ene =np.zeros(Nx)
    Pr = np.roll(pressure, 1)
    Pl = np.roll(pressure, -1)
    mom = U*velocity
    Er = np.roll(Energy, 1)
    El = np.roll(Energy, -1)
    for i in range(len(velocity)):
        if velocity[i] >=0:
            # Calculate flux at each interface
            F_rho[i] = 0.5 * velocity[i] * (Ur[i] + U[i]) - 0.5*np.abs(velocity[i]) * (Ur[i] - U[i])
            F_mom[i] = 0.5 * velocity[i] * (np.roll(F_rho, 1)[i] + F_rho[i] ) - 0.5*np.abs(velocity[i]) * (np.roll(F_rho, 1)[i] - F_rho[i]) + (Pr[i] - pressure[i])
            F_ene[i] = 0.5 * velocity[i] * (Er[i]+Pr[i] + Energy[i]+pressure[i]) - 0.5*np.abs(velocity[i]) * (Er[i]+Pr[i] - Energy[i]-pressure[i])
        else:
            # Calculate flux at each interface
            F_rho[i] = 0.5 * velocity[i] * (Ul[i] + U[i]) - 0.5*np.abs(velocity[i]) * (U[i] - Ul[i])
            F_mom[i] = 0.5 * velocity[i] * (np.roll(F_rho, 1)[i] + F_rho[i] ) - 0.5*np.abs(velocity[i]) * (F_rho[i] - np.roll(F_rho, -1)[i]) + (pressure[i] - Pl[i])
            F_ene[i] = 0.5 * velocity[i] * (El[i]+Pl[i] + Energy[i]+pressure[i]) - 0.5*np.abs(velocity[i]) * (Energy[i]+pressure[i] - El[i]-Pl[i])

    for i in range(len(velocity)):
        if velocity[i] >=0:
            # Calculate Residual
            Fr = np.roll(F_rho, 1)
            Fr_mom = np.roll(F_mom, 1)
            Fr_ene = np.roll(F_ene, 1)
            R[i] = F_rho[i] - Fr[i]
            R_mom[i] = F_mom[i] - Fr_mom[i]
            R_ene[i] = F_ene[i] - Fr_ene[i]

        else:
            # Calculate Residual
            Fl = np.roll(F_rho, -1)
            Fl_mom = np.roll(F_mom, -1)
            Fl_ene = np.roll(F_ene, -1)
            R[i] = Fl[i] - F_rho[i]
            R_mom[i] = Fl_mom[i] - F_mom[i]
            R_ene[i] = Fl_ene[i] - F_ene[i]

    U = U - (dt/dx)*R
    # Prevent 0 density
    for i in range(len(U)):
        if U[i] == 0.0:
            U[i] = 1e-6
    mom = mom - (dt/dx)*R_mom
    Energy = Energy - (dt/dx)*R_ene
    velocity = mom/U

    return U, velocity, Energy

def HLL(U, velocity,dt, dx, Nx, F_rho, F_mom, F_ene, pressure, Energy):
    Ur = np.roll(U, 1)
    Ul = np.roll(U, -1)
    R =np.zeros(Nx)
    R_mom =np.zeros(Nx)
    R_ene =np.zeros(Nx)
    Pr = np.roll(pressure, 1)
    Pl = np.roll(pressure, -1)
    mom = U*velocity
    Er = np.roll(Energy, 1)
    El = np.roll(Energy, -1)
    for i in range(len(velocity)):
        if velocity[i] >=0:
            # Calculate flux at each interface
            F_rho[i] = 0.5 * velocity[i] * (Ur[i] + U[i]) - 0.5*np.abs(velocity[i]) * (Ur[i] - U[i])
            F_mom[i] = 0.5 * velocity[i] * (np.roll(F_rho, 1)[i] + F_rho[i] ) - 0.5*np.abs(velocity[i]) * (np.roll(F_rho, 1)[i] - F_rho[i]) + (Pr[i] - pressure[i])
            F_ene[i] = 0.5 * velocity[i] * (Er[i]+Pr[i] + Energy[i]+pressure[i]) - 0.5*np.abs(velocity[i]) * (Er[i]+Pr[i] - Energy[i]-pressure[i])
        else:
            # Calculate flux at each interface
            F_rho[i] = 0.5 * velocity[i] * (Ul[i] + U[i]) - 0.5*np.abs(velocity[i]) * (U[i] - Ul[i])
            F_mom[i] = 0.5 * velocity[i] * (np.roll(F_rho, 1)[i] + F_rho[i] ) - 0.5*np.abs(velocity[i]) * (F_rho[i] - np.roll(F_rho, -1)[i]) + (pressure[i] - Pl[i])
            F_ene[i] = 0.5 * velocity[i] * (El[i]+Pl[i] + Energy[i]+pressure[i]) - 0.5*np.abs(velocity[i]) * (Energy[i]+pressure[i] - El[i]-Pl[i])

    for i in range(len(velocity)):
        if velocity[i] >=0:
            # Calculate Residual
            Fr = np.roll(F_rho, 1)
            Fr_mom = np.roll(F_mom, 1)
            Fr_ene = np.roll(F_ene, 1)
            R[i] = F_rho[i] - Fr[i]
            R_mom[i] = F_mom[i] - Fr_mom[i]
            R_ene[i] = F_ene[i] - Fr_ene[i]

        else:
            # Calculate Residual
            Fl = np.roll(F_rho, -1)
            Fl_mom = np.roll(F_mom, -1)
            Fl_ene = np.roll(F_ene, -1)
            R[i] = Fl[i] - F_rho[i]
            R_mom[i] = Fl_mom[i] - F_mom[i]
            R_ene[i] = Fl_ene[i] - F_ene[i]

    U = U - (dt/dx)*R
    # Prevent 0 density
    for i in range(len(U)):
        if U[i] == 0.0:
            U[i] = 1e-6
    mom = mom - (dt/dx)*R_mom
    Energy = Energy - (dt/dx)*R_ene
    velocity = mom/U

    return U, velocity, Energy


# Loop until t > tfinal
t = 0
n = 0
while (t < tfinal and n < 100) :

    # Calculate flux at each interface
    rho, velocity, Energy = upwind_flux(rho, velocity, dt, dx, Nx, F, F2, F3, pressure, Energy) # Forward Euler step

    # Derive primitives
    e = (Energy - 0.5*rho*velocity**2)/rho
    pressure = e*(gamma - 1)*rho

    # Time increment
    dt = CFL*dx/np.max(np.abs(velocity))
    t = t + dt
    n += 1

    #plt.ylim(-5, 5)
    #plt.plot(x[1:], velocity)
    plt.plot(x[1:], rho)
    plt.draw()
    plt.pause(0.0001)
    plt.clf()
