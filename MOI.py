from math import sin, cos, sqrt
import numpy as np
# import Shear Center as sh

""""
def f(theta):
    return -0.08*np.cos(theta)

def g(theta):
    return (-0.08+0.0011)*np.cos(theta)
"""

def locBooms():
                                                            # Given Parameters
    circ = 1.116
    Ca = 0.505  # m
    ha = 0.16  # m
    arc = circ/11
    le = arc/2
    theta = np.arctan((ha/2)/(Ca - (ha/2)))
    z2 = (ha/2) - (ha/2) * np.cos((arc/(ha/2)))
    y2 = (ha/2) * np.sin((arc/(ha/2)))
    zz = []
    yy = []
                                                            # z-coordinates
    zz.append(0)
    zz.append(z2)
    for i in range(0,4):
        zz.append(Ca - (le + 3*arc - i*arc)*np.cos(theta))
    for j in range(0,4):
        zz.append(Ca - (le + j * arc) * np.cos(theta))
    zz.append(z2)
                                                            # y-coordinates
    yy.append(0)
    yy.append(y2)
    for i in range(0, 4):
        yy.append((le + 3 * arc - i * arc) * np.sin(theta))
    for j in range(0, 4):
        yy.append(-(le + j * arc) * np.sin(theta))
    yy.append(y2)
                                                            # Merging of y and z coordinates in a list
    boomLoc = []
    for k in range(len(yy)):
        boomLoc.append([zz[k], yy[k]])

    # print(yy)
    # print(zz)
    print(boomLoc)
    # return yy
    # return zz
    return boomLoc



def centroid(boomLoc):                  ################## Circular Section Is Not Taken Into Account Yet ###################
    # Given Parameters
    y = 0
    skt = 0.0011
    spt = 0.0024
    Ca = 0.505
    ha = 0.16
    theta = np.arctan((ha / 2) / (Ca - (ha / 2)))
    l = sqrt((ha/2)**2 + (Ca - (ha/2))**2)

    boomArea = 3.456 * 10**(-5)    # [m^2]
    Az = []
    Azsk = ha * spt * (ha/2) + 2 * l * skt * (Ca - (l/2) * np.cos(theta))
    area = 11 * boomArea + ha * spt + 2 * l * skt

    for i in range(0,11):
        Az.append(boomArea * boomLoc[i][0])


    # location of the centroid is given by (z, y) measured from Leading Edge
    z = (sum(Az) + Azsk)/area
    print("z = ", z, "[m]")
    print("y = ", y, "[m]")
    return boomArea, z, y


def inertiaZZ(boomArea):                ################## Circular Section Is Not Taken Into Account Yet ###################
    d = []
    t = 0.0011
    beta = 0.186058177
    a = sqrt(0.08**2 + (0.505 - 0.08)**2)
    I_z_skin = 2 * ((t * a**3 * (sin(beta))**2)/12 + t * a * ((a/2) * sin(beta))**2)    # Still need to add the circular section

    I_zz = 0

    for i in range(0,11):                               # MOI of the booms
        d.append(boomLoc[i][1])
        I_zz = I_zz + boomArea * (d[i])**2

    I_zz = I_zz + I_z_skin                              # MOI_booms + MOI_skin
    print("I_zz =", I_zz, "[m^4]")
    return I_zz


def inertiaYY(boomArea, z):                         ################## Circular Section Is Not Taken Into Account Yet ###################
    d = []
    t = 0.0011
    beta = 0.186058177
    a = sqrt(0.08**2 + (0.505 - 0.08)**2)
    I_y_skin = 2 * ((t * a**3 * (cos(beta))**2)/12 + t * a * (0.505 - z - (a/2) * cos(beta))**2)   # Still need to add MOI of the circular section

    I_yy = 0

    for i in range(0, 11):                          # MOI of the booms
        d.append(z - boomLoc[i][1])
        I_yy = I_yy + boomArea * (d[i])**2

    I_yy = I_yy + I_y_skin                           # MOI_booms + MOI_skin
    print("I_yy =", I_yy, "[m^4]")
    return I_yy



boomLoc = locBooms()
centroid(boomLoc)
boomArea, z, y = centroid(boomLoc)
inertiaZZ(boomArea)
inertiaYY(boomArea, z)