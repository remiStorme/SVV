import Internal_Forces_and_Deflections as ifd
from MOI import Iz, Iy
import numpy as np
import math as m

'''
x_mz = ifd.x_mz
y_mz = ifd.y_mz

x_my = ifd.x_my
y_my = ifd.y_my
'''

def normStresses():
    x_mz = ifd.x_mz
    y_mz = ifd.y_mz

    x_my = ifd.x_my
    y_my = ifd.y_my

    sigmax = []
    #print(x_mz)
    #print(y_mz)
    #print(x_my)
    #print(y_my)

    for i in range(len(x_mz)):
        sigmax.append((y_my[i] * x_mz[i])/Iy + (y_mz[i] * x_my[i])/Iz)

    sigmax = [list(val) for val in sigmax]
    print(sigmax)

normStresses()