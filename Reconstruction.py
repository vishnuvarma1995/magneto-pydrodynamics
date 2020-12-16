# Methods to reconstruct interfaces to solve Riemann
# problem
import numpy as np

# HLLE Solver with choice of reconstruction
def HLLE_reconstruction_1D(q1,q2,q3,v,Cs, N, dx, reconstruction='lin',linearmethod='centered', slopelimiter=None):
    # Initialize arrayz
    F_rhoR = np.zeros(N)
    F_rhoL = np.zeros(N)
    F_momR = np.zeros(N)
    F_momL = np.zeros(N)
    F_eneR = np.zeros(N)
    F_eneL = np.zeros(N)

    # Calculate HLLE Riemann Solver shock/rarefaction velocities
    lambdaL = np.zeros(N)
    lambdaR = np.zeros(N)
    for i in range(0, N-1):
        lambdaL[i] = v[i] - Cs[i]
        lambdaR[i] = v[i+1] + Cs[i+1]

    # Piecewise constant reconstruction
    if reconstruction=='const':
        for i in range(1, N-1):
            if lambdaL[i] > 0.0:
                F_rhoR[i] = q1[i]*v[i]
                F_rhoL[i] = q1[i-1]*v[i]
                F_momR[i] = q2[i]*v[i]
                F_momL[i] = q2[i-1]*v[i]
                F_eneR[i] = q3[i]*v[i]
                F_eneL[i] = q3[i-1]*v[i]
            elif lambdaR[i] < 0.0:
                F_rhoR[i] = q1[i+1]*v[i]
                F_rhoL[i] = q1[i]*v[i]
                F_momR[i] = q2[i+1]*v[i]
                F_momL[i] = q2[i]*v[i]
                F_eneR[i] = q3[i+1]*v[i]
                F_eneL[i] = q3[i]*v[i]
            else:
                F_rhoR[i]  = (lambdaR[i]*q1[i]*v[i] - lambdaL[i]*q1[i+1]*v[i+1] + lambdaR[i]*lambdaL[i]*(q1[i+1]-q1[i]))/(lambdaR[i] - lambdaL[i])
                F_rhoL[i]  = (lambdaR[i]*q1[i-1]*v[i-1] - lambdaL[i]*q1[i]*v[i] + lambdaR[i]*lambdaL[i]*(q1[i]-q1[i-1]))/(lambdaR[i] - lambdaL[i])
                F_momR[i] = (lambdaR[i]*q2[i]*v[i] - lambdaL[i]*q2[i+1]*v[i+1] + lambdaR[i]*lambdaL[i]*(q2[i+1]-q2[i]))/(lambdaR[i] - lambdaL[i])
                F_momL[i] = (lambdaR[i]*q2[i-1]*v[i-1] - lambdaL[i]*q2[i]*v[i] + lambdaR[i]*lambdaL[i]*(q2[i]-q2[i-1]))/(lambdaR[i] - lambdaL[i])
                F_eneR[i] = (lambdaR[i]*q3[i]*v[i] - lambdaL[i]*q3[i+1]*v[i+1] + lambdaR[i]*lambdaL[i]*(q3[i+1]-q3[i]))/(lambdaR[i] - lambdaL[i])
                F_eneL[i] = (lambdaR[i]*q3[i-1]*v[i-1] - lambdaL[i]*q3[i]*v[i] + lambdaR[i]*lambdaL[i]*(q3[i]-q3[i-1]))/(lambdaR[i] - lambdaL[i])


    # Piecewise Linear reconstruction
    elif reconstruction=='lin':
        sigmaR1 = np.zeros(N)
        sigmaR2 = np.zeros(N)
        sigmaR3 = np.zeros(N)
        sigmaL1 = np.zeros(N)
        sigmaL2 = np.zeros(N)
        sigmaL3 = np.zeros(N)

        # Different choices of calculating sigma
        # Fromm's method
        if linearmethod == 'centered':
            for i in range(1, N-1):
                sigmaR1[i] = 0.5*(q1[i+1] - q1[i-1])/dx
                sigmaR2[i] = 0.5*(q2[i+1] - q2[i-1])/dx
                sigmaR3[i] = 0.5*(q3[i+1] - q3[i-1])/dx
                # sigmaR4[i] = 0.5*(p[i+1] - p[i-1])/dx
                sigmaL1[i] = 0.5*(q1[i+1] - q1[i-1])/dx
                sigmaL2[i] = 0.5*(q2[i+1] - q2[i-1])/dx
                sigmaL3[i] = 0.5*(q3[i+1] - q3[i-1])/dx
                # sigmaL4[i] = 0.5*(p[i+1] - p[i-1])/dx
        # Beam-Warming method
        elif linearmethod == 'upwind':
                for i in range(1, N-1):
                    sigmaR1[i] = (q1[i] - q1[i-1])/dx
                    sigmaR2[i] = (q2[i] - q2[i-1])/dx
                    sigmaR3[i] = (q3[i] - q3[i-1])/dx
                    sigmaL1[i] = (q1[i] - q1[i-1])/dx
                    sigmaL2[i] = (q2[i] - q2[i-1])/dx
                    sigmaL3[i] = (q3[i] - q3[i-1])/dx
        # Lax-Wendroff method
        elif linearmethod == 'downwind':
                for i in range(1, N-1):
                    sigmaR1[i] = (q1[i+1] - q1[i])/dx
                    sigmaR2[i] = (q2[i+1] - q2[i])/dx
                    sigmaR3[i] = (q3[i+1] - q3[i])/dx
                    sigmaL1[i] = (q1[i+1] - q1[i])/dx
                    sigmaL2[i] = (q2[i+1] - q2[i])/dx
                    sigmaL3[i] = (q3[i+1] - q3[i])/dx
        else:
            print('Other methods not implemented')

        # Calculate HLLE with linear reconstruction
        for i in range(1, N-1):
            if lambdaL[i] > 0.0:
                F_rhoR[i] = q1[i]*v[i] + 0.5*v[i]*sigmaR1[i]*(dx - v[i]*dx)
                F_rhoL[i] = q1[i-1]*v[i] + 0.5*v[i]*sigmaL1[i]*(dx - v[i]*dx)
                F_momR[i] = q2[i]*v[i] + 0.5*v[i]*sigmaR2[i]*(dx - v[i]*dx)
                F_momL[i] = q2[i-1]*v[i] + 0.5*v[i]*sigmaL2[i]*(dx - v[i]*dx)
                F_eneR[i] = q3[i]*v[i] + 0.5*v[i]*sigmaR3[i]*(dx - v[i]*dx)
                F_eneL[i] = q3[i-1]*v[i] + 0.5*v[i]*sigmaL3[i]*(dx - v[i]*dx)

            elif lambdaR[i] < 0.0:
                F_rhoR[i] = q1[i+1]*v[i] - 0.5*v[i]*sigmaR1[i]*(dx + v[i]*dx)
                F_rhoL[i] = q1[i]*v[i] - 0.5*v[i]*sigmaL1[i]*(dx + v[i]*dx)
                F_momR[i] = q2[i+1]*v[i] - 0.5*v[i]*sigmaR2[i]*(dx + v[i]*dx)
                F_momL[i] = q2[i]*v[i] - 0.5*v[i]*sigmaL2[i]*(dx + v[i]*dx)
                F_eneR[i] = q3[i+1]*v[i] - 0.5*v[i]*sigmaR3[i]*(dx + v[i]*dx)
                F_eneL[i] = q3[i]*v[i] - 0.5*v[i]*sigmaL3[i]*(dx + v[i]*dx)

            else:
                F_rhoR[i]  = (lambdaR[i]*q1[i]*v[i] - lambdaL[i]*q1[i+1]*v[i+1] + lambdaR[i]*lambdaL[i]*(q1[i+1]-q1[i]))/(lambdaR[i] - lambdaL[i])
                F_rhoL[i]  = (lambdaR[i]*q1[i-1]*v[i-1] - lambdaL[i]*q1[i]*v[i] + lambdaR[i]*lambdaL[i]*(q1[i]-q1[i-1]))/(lambdaR[i] - lambdaL[i])
                F_momR[i] = (lambdaR[i]*q2[i]*v[i] - lambdaL[i]*q2[i+1]*v[i+1] + lambdaR[i]*lambdaL[i]*(q2[i+1]-q2[i]))/(lambdaR[i] - lambdaL[i])
                F_momL[i] = (lambdaR[i]*q2[i-1]*v[i-1] - lambdaL[i]*q2[i]*v[i] + lambdaR[i]*lambdaL[i]*(q2[i]-q2[i-1]))/(lambdaR[i] - lambdaL[i])
                F_eneR[i] = (lambdaR[i]*q3[i]*v[i] - lambdaL[i]*q3[i+1]*v[i+1] + lambdaR[i]*lambdaL[i]*(q3[i+1]-q3[i]))/(lambdaR[i] - lambdaL[i])
                F_eneL[i] = (lambdaR[i]*q3[i-1]*v[i-1] - lambdaL[i]*q3[i]*v[i] + lambdaR[i]*lambdaL[i]*(q3[i]-q3[i-1]))/(lambdaR[i] - lambdaL[i])

    return F_rhoR, F_rhoL, F_momR, F_momL, F_eneR, F_eneL
