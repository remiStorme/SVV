from math import pi
from math import cos
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
plt.style.use('ggplot')


Nz = 81
Nx = 41
ca = 0.505
la = 1.611


def Coordinates():

    def thetaZ(Nz):
        thz = []
        for i in range(1,Nz+2):
            thz.append((i-1)*pi/Nz)
        return thz

    def ZCoordinate(ca,Nz,thz):
        z = []
        for i in range(0,Nz):
            z.append(-0.25*ca*((1-cos(thz[i]))+(1-cos(thz[i+1]))))
        return z

    def thetaX(Nx):
        thx = []
        for i in range(1,Nx+2):
            thx.append((i-1)*pi/Nx)
        return thx

    def XCoordinate(la,Nx,thx):
        x = []
        for i in range(0,Nx):
            x.append(0.25*la*((1-cos(thx[i]))+(1-cos(thx[i+1]))))
        return x

    ThetaZ = thetaZ(Nz)
    ZCoord = ZCoordinate(ca, Nz, ThetaZ)
    ThetaX = thetaZ(Nx)
    XCoord = XCoordinate(la, Nx, ThetaX)

    return ZCoord,XCoord


def FileReader(path):

    rows = []
    with open(path,"r+") as f:
        lines = f.readlines()
        i = 0
        for line in lines:
            line = line.strip()
            line = line.split(",")
            line = [float(val) for val in line]
            rows.append(line)

    mat = np.array(rows)
    return mat


def plotter(XCoord,mat,idx):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Chordwise Location')
    ax.set_ylabel('Value')
    ax.set_title("Aerodynamic Data for Chord Slice")
    #x = np.arange(mat.shape[1])
    #plt.plot(XCoord,mat[idx,:],"r+")
    #plt.show()

mat = FileReader("C:/Users/Paul Simon Schön/Downloads/aerodynamicloadf100.dat")

X,Z = Coordinates()
plotter(XCoord,mat,0)

''''


def surface_plot (XCoord,ZCoord,matrix, **kwargs):
    # acquire the cartesian coordinate matrices from the matrix
    # x is cols, y is rows
    (x, y) = np.meshgrid(81, 41)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(XCoord, ZCoord, matrix, **kwargs)
    return (fig, ax, surf)

print(mat.shape)
(fig2, ax2, surf) = surface_plot(ZCoord,XCoord,mat,cmap=plt.cm.coolwarm)

fig2.colorbar(surf)

ax2.set_xlabel('Z (cols)')
ax2.set_ylabel('X (rows)')
ax2.set_zlabel('values')

plt.show()

''''