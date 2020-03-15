import Internal_Forces_and_Deflections as ifd
from MOI import Iz, Iy, z
import matplotlib.pyplot as plt
import numpy as np
import math as m

'''
x_mz = ifd.x_mz
y_mz = ifd.y_mz

x_my = ifd.x_my
y_my = ifd.y_my
'''
Ca = 0.505
ha = 0.161
zz = np.linspace(0, Ca, 40)
yy = np.linspace(0, ha, 40)

def normStresses():
    x_mz = ifd.x_mz     #np.array(np.linspace(0, la, 40))
    y_mz = ifd.y_mz

    x_my = ifd.x_my     #np.array(np.linspace(0, la, 40))
    y_my = ifd.y_my

    sigmax = []
    #print(x_mz)
    #print(y_mz)
    #print(x_my)
    #print(y_my)

    for i in range(len(zz)):
        sigmax.append((y_my[i] * (zz[i] - z))/Iy + (y_mz[i] * yy)/Iz)

    sigmax = [list(val) for val in sigmax]
    #print(sigmax)

    plt.plot(x_mz, sigmax)
    plt.xlabel('x-position (m)')
    plt.ylabel('Sigma x (Pa)')
    plt.show()
    print(sigmax)
    return sigmax
normStresses()