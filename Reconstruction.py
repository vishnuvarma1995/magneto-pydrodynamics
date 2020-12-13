# Methods to reconstruct interfaces to solve Riemann
# problem

def piecewise_constant():
    pass

def piecewise_linear():
    pass

# Old 1D
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


    # Piecewise Linear reconstruction
    sigmaR1 = np.zeros(N)
    sigmaR2 = np.zeros(N)
    sigmaR3 = np.zeros(N)
    sigmaR4 = np.zeros(N)
    sigmaL1 = np.zeros(N)
    sigmaL2 = np.zeros(N)
    sigmaL3 = np.zeros(N)
    sigmaL4 = np.zeros(N)
    # Fromm's method
    for i in range(1, N-1):
        sigmaR1[i] = 0.5*(q1[i+1] - q1[i-1])/dx
        sigmaR2[i] = 0.5*(q2[i+1] - q2[i-1])/dx
        sigmaR3[i] = 0.5*(q3[i+1] - q3[i-1])/dx
        sigmaR4[i] = 0.5*(p[i+1] - p[i-1])/dx
        sigmaL1[i] = 0.5*(q1[i+1] - q1[i-1])/dx
        sigmaL2[i] = 0.5*(q2[i+1] - q2[i-1])/dx
        sigmaL3[i] = 0.5*(q3[i+1] - q3[i-1])/dx
        sigmaL4[i] = 0.5*(p[i+1] - p[i-1])/dx
    for i in range(1, N-1):
        if v[i] > 0.0:
            F_rhoR[i] = q1[i]*v[i] + 0.5*v[i]*sigmaR1[i]*(dx - v[i]*dx)
            F_rhoL[i] = q1[i-1]*v[i] + 0.5*v[i]*sigmaL1[i]*(dx - v[i]*dx)
            F_mom1R[i] = q2[i]*v[i] + 0.5*v[i]*sigmaR2[i]*(dx - v[i]*dx)
            F_mom1L[i] = q2[i-1]*v[i] + 0.5*v[i]*sigmaL2[i]*(dx - v[i]*dx)
            F_ene1R[i] = q3[i]*v[i] + 0.5*v[i]*sigmaR3[i]*(dx - v[i]*dx)
            F_ene1L[i] = q3[i-1]*v[i] + 0.5*v[i]*sigmaL3[i]*(dx - v[i]*dx)
            F_ene2R[i] = p[i]*v[i] #+ 0.5*v[i]*sigmaR4*(dx - v[i]*dx)
            F_ene2L[i] = p[i-1]*v[i] #+ 0.5*v[i]*sigmaL4*(dx - v[i]*dx)
        elif v[i] < 0.0:
            F_rhoR[i] = q1[i+1]*v[i] - 0.5*v[i]*sigmaR1[i]*(dx + v[i]*dx)
            F_rhoL[i] = q1[i]*v[i] - 0.5*v[i]*sigmaL1[i]*(dx + v[i]*dx)
            F_mom1R[i] = q2[i+1]*v[i] - 0.5*v[i]*sigmaR2[i]*(dx + v[i]*dx)
            F_mom1L[i] = q2[i]*v[i] - 0.5*v[i]*sigmaL2[i]*(dx + v[i]*dx)
            F_ene1R[i] = q3[i+1]*v[i] - 0.5*v[i]*sigmaR3[i]*(dx + v[i]*dx)
            F_ene1L[i] = q3[i]*v[i] - 0.5*v[i]*sigmaL3[i]*(dx + v[i]*dx)
            F_ene2R[i] = p[i+1]*v[i] #- 0.5*v[i]*sigmaR4[i]*(dx + v[i]*dx)
            F_ene2L[i] = p[i]*v[i] #- 0.5*v[i]*sigmaL4[i]*(dx + v[i]*dx)
        else:
            F_rhoR[i]  = 0.0
            F_rhoL[i]  = 0.0
            F_mom1R[i] = 0.0
            F_mom1L[i] = 0.0
            F_ene1R[i] = 0.0
            F_ene1L[i] = 0.0
            F_ene2R[i] = 0.0
            F_ene2L[i] = 0.0
