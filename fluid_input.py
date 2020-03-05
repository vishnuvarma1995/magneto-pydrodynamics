# File structure based on THE_ARGO
# by Vishnu Varma

class TheFluid (object):

    def __init__ (self, Nx, Ny, Nz, Dx, Dy, Dz):
        # Number of zones along x, y and z 
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz

        # Width of cell in x, y and z
        self.Dx = Dx
        self.Dy = Dy
        self.Dz = Dz

# end class TheFluid