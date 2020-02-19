

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
        Az.append(boomLoc[i][0])

    z = sum(Az)/area
    print(z)


I_zz=1  #temp value for lorenzas code