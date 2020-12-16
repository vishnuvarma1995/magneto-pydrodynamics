import EulerAdvection as hydro_1d


'''
Date: 16/12/2020
Setup model initial conditions, dimensionality and numerical methods to use
'''

DIMN = 1 # Dimensionality of problem (1 or 2 dimensions)
Nx = 500 # Grid points in x-direction
Ny = 0 # Grid points in y -direction
CFL = 0.7
tfinal = 0.25
gamma = 1.66 # For gamma law pressure, 1.66 for ideal gas

# Initialize arrays
rho = np.zeros(Nx)
energy = np.zeros(Nx)
pressure = np.zeros(Nx)
velocity = np.zeros(Nx)

# Test problems
# Sod's shock tube problem
pressure[:int(Nx/2)] = 1.0
pressure[int(Nx/2):] = 0.1
rho[:int(Nx/2)] = 1.0
rho[int(Nx/2):] = 0.125
velocity[:int(Nx/2)] = 0.0
velocity[int(Nx/2):] = 0.0

#Blast Wave Problem
rho[:] = 1.0
velocity[:] = 0.0
pressure[:int(Nx/10)] = 1000
pressure[int(Nx/10):int(9*Nx/10)] = 0.1
pressure[int(9*Nx/10):] = 100

# 123 Problem
rho[:] = 1.0
velocity[:int(Nx/2)] = -2.0
velocity[int(Nx/2):] = 2.0
pressure[:] = 0.4
