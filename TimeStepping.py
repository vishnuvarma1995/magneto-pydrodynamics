# Choices of methods for performing numberical timestepping

def upwind(u, fR, fL, dt, dx):
    return u - (dt/dx)*(fR - fL)


def LaxFriedrichs(uR, uL, fR, fL, dt, dx):
    return 0.5*(uR + uL) - 0.5*(dt/dx)*(fR - fL)
