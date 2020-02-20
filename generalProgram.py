import MOI
import AEload
import SC_copy

boomLoc = locBooms()
boomArea, z, y = centroid(boomLoc)
mat = FileReader("C:/Users/Paul Simon Schön/Downloads/aerodynamicloadf100.dat")
X, Z = Coordinates()



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






