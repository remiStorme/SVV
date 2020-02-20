import MOI
import AEload.py

boomLoc = locBooms()
boomArea, z, y = centroid(boomLoc)
inertiaZZ(boomArea)
inertiaYY(boomArea, z)



class aileron:
    def __init__(self,moi_zz,moi_xx,shear_c,ca,la):
        self.moi_zz = moi_zz
        self.shear_c = shear_c
        self.centroid = centroid
        self.moi_xx = moi_xx
        self.ca = 0.505 #m
        self.la = 1.611 #m
        self.x_pos = []
        self.z_pos = []
        self.hinge_1 = 0.125 #m
        self.hinge_2 = 0.498 #m
        self.hinge_3 = 1.494 #m
        self.xa = 0.245 #m
        self.P = 49.2 #kN



