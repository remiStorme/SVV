#import MOI
from AEload import FileReader,Coordinates
#import SC
from interpolator_Cubic import Interpolate

#under construction

#boomLoc = locBooms()
#boomArea, z, y = centroid(boomLoc)
mat = FileReader("C:/Users/Paul Simon Schön/Downloads/aerodynamicloadf100.dat")
X, Z = Coordinates()


class chord_slice:

    def __init__(self,distr,zpos):
        self.z_pos = zpos
        self.distr = distr

        def getFunc(self):
            c = Interpolate(self.distr,self.z_pos)
            func, integral = c.spline_natural()
            return func, integral

        func, integral = self.getFunc()
        self.func = func
        self.integral = integral

    def shearForceLoc(self):
        idx = 0
        maxval = max(self.func)
        for i in range(len(self.z_pos)):
            if self.func[i] != maxval:
                continue
            else:
                idx = i
        else:
            print("Error")
            return
        return self.z_pos[i],idx

    def aero_Func(self):
        pass

    def Torque(self):
        pass

    def shearFlow(self):
        pass

    def NormalStress(self):
        pass





class aileron:

    def __init__(self,moi_zz,moi_xx,shear_c,ca,la):
        #self.moi_zz = inertiaZZ(boomArea)
        #self.shear_c = shear_c
        #self.centroid_z = z
        #self.centroid_y = y
        #self.moi_yy = inertiaYY(boomArea, z)
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

    def genSlices(self):
        for i in range(len(self.mat)):
            self.slices.append(chord_slice(m[i,:],self.Zcoord))




    def Moment_Distribution(self):
        pass

    def VanMieses(self):
        pass

    def deflection(self):
        pass

    def twist(self):
        pass










test = aileron(0,0,0,0,0)
test_list = test.genSlices()
# print(test_list)


