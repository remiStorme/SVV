

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

    for i in range(0,11):
        Az = []
        Az.append(boomArea * boomLoc[i][0])

    # location of the centroid is given by (z, y) measured from Leading Edge
    z = sum(Az)/area
    print("z = ", z)
    print("y = ", y)
    return boomLoc, boomArea, z, y


def inertiaZZ(boomArea):
    d = []
    I_zz = 0

    for i in range(0,11):
        d.append(boomLoc[i][1])
        I_zz = I_zz + boomArea * (d[i])**2
    print("I_zz is", I_zz)
    return I_zz


def inertiaYY(boomArea, z):
    d = []
    I_yy = 0

    for i in range(0, 11):
        d.append(z - boomLoc[i][1])
        I_yy = I_yy + boomArea * (d[i])**2
    print("I_yy is", I_yy)
    return I_yy

centroid()
boomLoc, boomArea, z, y = centroid()
inertiaZZ(boomArea)
inertiaYY(boomArea, z)
