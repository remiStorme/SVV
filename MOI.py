from math import sin, cos, sqrt
import numpy as np

class CrossSection:

    def __init__(self):
        self.circ = 1.116
        self.Ca = 0.505  # m
        self.ha = 0.16  # m
        self.arc = self.circ / 11
        self.le = self.arc / 2
        self.theta = np.arctan((self.ha / 2) / (self.Ca - (self.ha / 2)))
        self.skt = 0.0011
        self.spt = 0.0024
        self.l = sqrt((self.ha / 2) ** 2 + (self.Ca - (self.ha / 2)) ** 2)
        self.boomArea = 3.456 * 10 ** (-5)

    def locBooms(self):
                                                                # Given Parameter
        z2 = (self.ha/2) - (self.ha/2) * np.cos((self.arc/(self.ha/2)))
        y2 = (self.ha/2) * np.sin((self.arc/(self.ha/2)))
        zz = []
        yy = []
                                                                # z-coordinates
        zz.append(0)
        zz.append(z2)
        for i in range(0,4):
            zz.append(self.Ca - (self.le + 3*self.arc - i*self.arc)*np.cos(self.theta))
        for j in range(0,4):
            zz.append(self.Ca - (self.le + j * self.arc) * np.cos(self.theta))
        zz.append(z2)
                                                                # y-coordinates
        yy.append(0)
        yy.append(y2)
        for i in range(0, 4):
            yy.append((self.le + 3 * self.arc - i * self.arc) * np.sin(self.theta))
        for j in range(0, 4):
            yy.append(-(self.le + j * self.arc) * np.sin(self.theta))
        yy.append(-y2)
                                                                # Merging of y and z coordinates in a list
        boomLoc = []
        for k in range(len(yy)):
            boomLoc.append([zz[k], yy[k]])

        # print(yy)
        # print(zz)
        # print(boomLoc)
        # return yy
        # return zz
        return boomLoc



    def centroid(self):
        # Given Parameters
        y = 0
        boomLoc = self.locBooms()
        Az = []
        Azsk = self.ha * self.spt * (self.ha/2) + 2 * self.l * self.skt * (self.Ca - (self.l/2) * np.cos(self.theta)) + ((self.ha/2) - (2 * (self.ha/2))/np.pi) * np.pi * (self.ha/2) * self.skt
        area = 11 * self.boomArea + self.ha * self.spt + 2 * self.l * self.skt + np.pi * (self.ha/2) * self.skt

        for i in range(0,11):
            Az.append(self.boomArea * boomLoc[i][0])

        # location of the centroid is given by (z, y) measured from Leading Edge
        z = (sum(Az) + Azsk)/area
        # print("The coordinates of the centroid measured from the Leading Edge are (", z, ",", y, ") [m]")
        #print("Centroid is:", z,y)
        return z, y


    def inertiaZZ(self):
        d = []
        boomLoc = self.locBooms()
        I_z_circ = (np.pi * self.skt * (self.ha / 2)**3) / 2
        I_z_skin = 2 * ((self.skt * self.l**3 * (sin(self.theta))**2)/12 + self.skt * self.l * ((self.l/2) *
        sin(self.theta))**2) + (self.spt * self.ha**3)/12 + I_z_circ  # Still need to add the circular section

        I_zz = 0

        for i in range(0,11):                               # MOI of the booms
            d.append(boomLoc[i][1])
            I_zz = I_zz + self.boomArea * (d[i])**2


        I_zz = I_zz + I_z_skin                              # MOI_booms + MOI_skin
        # print("I_zz =", I_zz, "[m^4]")
        return I_zz


    def inertiaYY(self):
        d = []
        z,_ = self.centroid()
        boomLoc = self.locBooms()
        I_y_circ = (np.pi * self.skt * (self.ha / 2) ** 3) / 2 - np.pi * (self.ha / 2) * self.skt * ((2 * (self.ha / 2))
        / np.pi) ** 2 + np.pi * (self.ha / 2) * self.skt * (z - ((self.ha/2) - (2 * (self.ha/2))/np.pi))**2
        I_y_skin = 2 * ((self.skt * self.l**3 * (cos(self.theta))**2)/12 + self.skt * self.l * (0.505 - z - (self.l/2) *
        cos(self.theta))**2) + (self.spt**3 * self.ha)/12 + self.ha * self.spt * (z - self.ha)**2 + I_y_circ

        I_yy = 0

        for i in range(0, 11):                          # MOI of the booms
            d.append(z - boomLoc[i][1])
            I_yy = I_yy + self.boomArea * (d[i])**2

        I_yy = I_yy + I_y_skin                          # MOI_booms + MOI_skin
        # print("I_yy =", I_yy, "[m^4]")
        return I_yy


                                                        # Calling the functions
cs = CrossSection()
locbooms = cs.locBooms()
centroid = cs.centroid()
Iz = cs.inertiaZZ()
Iy = cs.inertiaYY()

# The error of the values of Izz and Iyy compared to the ones from the verification model

e_yy = ((Iy - 4.5943507864451845 * 10**(-5))/(4.5943507864451845 * 10**(-5))) * 100     # %
e_zz = ((4.753851442684436 * 10**(-6) - Iz)/(4.753851442684436 * 10**(-6))) * 100       # %

print('Iyy =', Iy,'and Izz =', Iz)
print('The error on Iyy is', e_yy,'%')
print('The error on Izz is', e_zz,'%')
z,y = centroid

e_z = ((z - 0.20362591085157106)/(0.20362591085157106)) * 100           # %
print('The error on z is', z,'%')

'''
plt.figure()
plt.scatter(zcoord,y, marker="o", c= tau_yz)
plt.title("Shear stress in the cross section at ")
plt.xlabel("z")
plt.ylabel("y")
plt.colorbar()
plt.axis('equal')
plt.show
'''