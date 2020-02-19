from math import sin, cos, sqrt

def centroid():
    # Boom Coordinates measured from the Leading Edge in [m]
    zy1 = [0, 0]
    zy2 = [0.056158645, 0.076364847]
    zy3 = [0.156037587, 0.065687042]
    zy4 = [0.255741134, 0.046919316]
    zy5 = [0.35544468, 0.02815159]
    zy6 = [0.455148226, 0.009383863]
    zy7 = [0.455148226, -0.009383863]
    zy8 = [0.35544468, -0.02815159]
    zy9 = [0.255741134, -0.046919316]
    zy10 = [0.156037587, -0.065687042]
    zy11 = [0.056158645, -0.076364847]

    # Centroid Coordinates
    z = 0
    y = 0

    boomLoc = [zy1, zy2, zy3, zy4, zy5, zy6, zy7, zy8, zy9, zy10, zy11]
    boomArea = 3.456 * 10**(-5)    # [m^2]
    area = 11 * boomArea
    Az = []

    for i in range(0,11):
        Az.append(boomArea * boomLoc[i][0])


    # location of the centroid is given by (z, y) measured from Leading Edge
    z = sum(Az)/area
    # print("z = ", z, "[m]")
    # print("y = ", y, "[m]")
    return boomLoc, boomArea, z, y


def inertiaZZ(boomArea):
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
    # print("I_zz is", I_zz, "[m^4]")
    return I_zz


def inertiaYY(boomArea, z):
    d = []
    t = 0.0011
    beta = 0.186058177
    a = sqrt(0.08**2 + (0.505 - 0.08)**2)
    I_y_skin = 2 * ((t * a**3 * (cos(beta))**2)/12 + t * a * (0.505 - z - (a/2) * cos(beta))**2)       # Still need to add MOI of the circular section

    I_yy = 0

    for i in range(0, 11):                          # MOI of the booms
        d.append(z - boomLoc[i][1])
        I_yy = I_yy + boomArea * (d[i])**2

    I_yy = I_yy + I_y_skin                           # MOI_booms + MOI_skin
    # print("I_yy is", I_yy, "[m^4]")
    return I_yy

boomLoc, boomArea, z, y = centroid()
inertiaZZ(boomArea)
inertiaYY(boomArea, z)
