import numpy as np

# Options for approximate Riemann solvers

def HLLE(uL, uR, fL, fR, S_L, S_R, N):
    # Conserved Variables
    F  = np.zeros(N)
    U  = np.zeros(N)

    for i in range(1,N-2):
        if min(S_L[i], 0.0) > 0.0:
            F[i] = fL[i]
            U[i] = uL[i]

        elif max(S_R[i], 0.0) < 0.0:
            F[i] = fR[i]
            U[i] = uR[i]

        elif min(S_L[i], 0.0) <= 0.0 and max(S_R[i], 0.0) >= 0.0:
            F[i] = ( S_R[i]*fL[i] - S_L[i]*fR[i] + S_R[i]*S_L[i]*(uL[i] - uR[i]) )/(S_R[i] - S_L[i] + 1e-11)
            U[i] = ( S_R[i]*uR[i] - S_L[i]*uL[i] + fR[i] + fL[i] )/(S_R[i] - S_L[i] + 1e-11)
    # Compute L/R fluxes along the lines bm/bp: F_{L}-S_{L}U_{L}; F_{R}-S_{R}U_{R}

    # Add left boundary
    F[0]  = F[1]
    F[-2] = F[-3]
    F[-1] = F[-2]
    U[0]  = U[1]
    U[-2] = U[-3]
    U[-1] = U[-2]

    return U, F


def HLLE_ll(rL, rR, vL, vR, pL, pR, gamma, N):
    # Left State
    uL_1 = rL # density
    uL_2 = rL * vL # momentum
    uL_3 = pL/((gamma-1)) + 0.5*(uL_2**2)/uL_1 # Energy
    Cs_L = np.sqrt(np.abs(gamma*pL/rL))#Sound speed
    lambda_Lr = vL + Cs_L
    lambda_Ll = vL - Cs_L

    # Left Fluxes
    F = uL_2
    F2 = uL_2*vL + pL
    F3 = (uL_3 + pL)*vL

    F_L  = np.zeros(N)
    F2_L = np.zeros(N)
    F3_L = np.zeros(N)

    for i in range(1,N-2):
        if lambda_Ll[i] >= 0.0:
            F_L[i]  = F[i]
            F2_L[i] = F2[i]
            F3_L[i] = F3[i]

        elif lambda_Lr[i] <= 0.0:
            F_L[i]  = F[i+1]
            F2_L[i] = F2[i+1]
            F3_L[i] = F3[i+1]

        elif lambda_Ll[i] < 0.0 and lambda_Lr[i] > 0.0:
            F_L[i] = (lambda_Lr[i] *F[i] - lambda_Ll[i] *F[i+1] + lambda_Lr[i] *lambda_Ll[i] *(uL_1[i+1]-uL_1[i]))/(lambda_Lr[i]  - lambda_Ll[i] )
            F2_L[i] =(lambda_Lr[i] *F[i] - lambda_Ll[i] *F[i+1] + lambda_Lr[i] *lambda_Ll[i] *(uL_2[i+1]-uL_2[i]))/(lambda_Lr[i]  - lambda_Ll[i] )
            F3_L[i] =(lambda_Lr[i] *F[i] - lambda_Ll[i] *F[i+1] + lambda_Lr[i] *lambda_Ll[i] *(uL_3[i+1]-uL_3[i]))/(lambda_Lr[i]  - lambda_Ll[i] )

    # Right state
    uR_1 = rR # density
    uR_2 = rR * vR # momentum
    uR_3 = pR/((gamma-1)) + 0.5*(uR_2**2)/uR_1 # Energy
    Cs_R = np.sqrt(np.abs(gamma*pR/rR))#Sound speed
    lambda_Rr = vR + Cs_R
    lambda_Rl = vR - Cs_R
    # Right Fluxes
    F = uR_2
    F2 = uR_2*vR + pR
    F3 = (uR_3 + pR)*vR

    F_R  = np.zeros(N)
    F2_R = np.zeros(N)
    F3_R = np.zeros(N)

    for i in range(1,N-2):
        if lambda_Rl[i] >= 0.0:
            F_R[i]  = F[i]
            F2_R[i] = F2[i]
            F3_R[i] = F3[i]

        elif lambda_Rr[i] <= 0.0:
            F_R[i]  = F[i+1]
            F2_R[i] = F2[i+1]
            F3_R[i] = F3[i+1]

        elif lambda_Rl[i] < 0.0 and lambda_Rr[i] > 0.0:
            F_R[i] = (lambda_Rr[i] *F[i] - lambda_Rl[i] *F[i+1] + lambda_Rr[i] *lambda_Rl[i] *(uR_1[i+1]-uR_1[i]))/(lambda_Rr[i]  - lambda_Rl[i] )
            F2_R[i] =(lambda_Rr[i] *F[i] - lambda_Rl[i] *F[i+1] + lambda_Rr[i] *lambda_Rl[i] *(uR_2[i+1]-uR_2[i]))/(lambda_Rr[i]  - lambda_Rl[i] )
            F3_R[i] =(lambda_Rr[i] *F[i] - lambda_Rl[i] *F[i+1] + lambda_Rr[i] *lambda_Rl[i] *(uR_3[i+1]-uR_3[i]))/(lambda_Rr[i]  - lambda_Rl[i] )

    # % Wave speed estimates
    # RT = sqrt(rR/rL); % r = RT*rL;
    # u = (uL+RT*uR)/(1+RT);
    # SLm = min([ vnL-aL, vn-a, 0]);
    # SRp = max([ vnR+aR, vn+a, 0]);
    #
    # % Left and Right fluxes
    # FL=[rL*vnL; rL*vnL*uL + pL*nx; rL*vnL*vL + pL*ny; rL*vnL*HL];
    # FR=[rR*vnR; rR*vnR*uR + pR*nx; rR*vnR*vR + pR*ny; rR*vnR*HR];

    # Add left boundary
    F_R[0] = F_R[1]
    F2_R[0] = F2_R[1]
    F3_R[0] = F3_R[1]
    F_R[N-1] = F_R[N-2]
    F2_R[N-1] = F2_R[N-2]
    F3_R[N-1] = F3_R[N-2]

    F_L[0] = F_L[1]
    F2_L[0] = F2_L[1]
    F3_L[0] = F3_L[1]
    F_L[N-1] = F_L[N-2]
    F2_L[N-1] = F2_L[N-2]
    F3_L[N-1] = F3_L[N-2]
    return F_R, F2_R, F3_R, F_L, F2_L, F3_L
