import MOI
import AEload
import SC_copy

#under construction

boomLoc = locBooms()
boomArea, z, y = centroid(boomLoc)
mat = FileReader("C:/Users/Paul Simon Schön/Downloads/aerodynamicloadf100.dat")
X, Z = Coordinates()


class chord_slice:

    def __init__(self,xpos,Z):
        self.x_pos = xpos
        self.zCoord = Z
        self.distr = []

    def shearForceLoc(self):

    def Max_aero_Force(self):

    def Mag_aero_Force(self):

    def aero_Func(self):

    def Torque(self):

    def shearFlow(self):

    def NormalStress(self):








chord_slice(self):

class aileron:

    def __init__(self,moi_zz,moi_xx,shear_c,ca,la):
        self.moi_zz = inertiaZZ(boomArea)
        self.shear_c = shear_c
        self.centroid_z = z
        self.centroid_y = y
        self.moi_yy = inertiaYY(boomArea, z)
        self.ca = 0.505 #m
        self.la = 1.611 #m
        self.x_pos = []
        self.z_pos = []
        self.hinge_1 = 0.125 #m
        self.hinge_2 = 0.498 #m
        self.hinge_3 = 1.494 #m
        self.xa = 0.245 #m
        self.P = 49.2 #kN
        self.mat = mat
        self.Xcoord = X
        self.ZCoord = Z
        self.slices = []

    def generateSlices(self):
        self.xlocs = []
        self.distributions = [mat[]]


    def Moment_Distribution(self):

    def VanMieses(self):

    def deflection(self):

    def twist(self):













