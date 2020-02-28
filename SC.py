######################## IMPORTS ########################################################

from numpy import *
from math import pi
from MOI import CrossSection
import Reaction_Forces_Aero_Load as rr
import AEload as ae
import Internal_Forces_and_Deflections as ifd


#HINGE AND ACTUATOR LOCATIONS: TEMPORARY VARIABLES
# x1 = 0.125  # m
# x2 = 0.498  # m
# x3 = 1.494  # m
# xaj     = x2 - xa/2
# xp      = x2 + xa/2



x_T,  y_T = ifd.x_T,  ifd.y_T #total torque on x axis around shear centre generated by ae load and reaction forces
x_sy, y_sy = ifd.x_sy, ifd.y_sy #x_sy is the slice at which sy is being evaluated and y_sy is the magnitude of Sy at
# that slice
x_sz, y_sz = ifd.x_sz, ifd.y_sz #x_sz is the slice at which sz is being evaluated and y_sz is the magnitude of Sz at
# that slice
lift_ae = ae.w_x #resultant lift per cross section

#vector x = [Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Ar, C1, C2, C3, C4, C5]
uks = rr.uks


######################## GEOMETRY F100 ########################################################

Ca = 0.505  # m
la = 1.611  # m
x1 = 0.125  # m
x2 = 0.498  # m
x3 = 1.494  # m
xa = 0.25  # m
ha = 0.16  # m
tsk = 1.1 / 1000  # m
tsp = 2.4 / 1000  # m
tst = 1.2 / 1000  # m
hst = 13. / 1000  # m
wst = 17. / 1000  # m
nst = 11  # m
d1 = 0.00389  # m
d3 = 0.01245  # m
theta = radians(30)  # rad
P = 49.2 * 1000  # N
r = ha / 2  # m
ladiag = sqrt((Ca - r) ** 2 + (r) ** 2)
alpha = arctan(r/(Ca-r))
cs = CrossSection()
boomLoc = cs.locBooms()
boomArea = 3.456 * 10 ** (-5)
z, y = cs.centroid()
I_zz = 4.686*10**-6 #from verification
I_yy = cs.inertiaZZ()

######################## DISCRETIZATION PARAMETERS ################################################

N = 100 # !!!!!!!! do evenly discretized cross section, not like it is now
nbooms = 11 #!!!!!!!!! include boom 0 for when Sx will be introduced

######################## INTEGRATION TOOLS #####################################################

def getszyf():
    s = []
    s0 = linspace(0, (pi / 2),N)  # (start,stop,number of steps). s is the angle going from 0 to boom1, so its actually
    # theta and the integral infact is in dtheta not ds
    s1 = linspace(0, r, N)
    s2 = linspace(0, ladiag, N)
    s3 = s2
    s4 = linspace(0, -r, N)
    s5 = s0
    s.extend((s0, s1, s2, s3, s4, s5))

    y = []
    y0 = []
    for i in range(len(s[0])):
        y0.append(r * sin(s[0][i]))
    y1 = s[1]
    y2 = []
    for i in range(len(s[2])):
        y2_i = ((ladiag - s[2][i]) / ladiag) * r
        y2.append(y2_i)
    y3 = [element - r for element in y2]
    y4 = [element * -1 for element in y1]
    y5=  []
    for i in range(len(y0)):
        y5.append(-y0[N - 1 - i])
    y.extend((y0, y1, y2, y3, y4, y5))

    z= []
    z0=[]
    for i in range(len(s[0])):
        z0.append(-(r-cos(s[0][i])*r))
    z1= [r]*N
    z2=[]
    for i in range(len(s[2])):
        z2.append(-(r+s[2][i]*cos(alpha)))
    z3=[]
    for i in range(len(s[3])):
        z3.append(-(Ca-s[3][i]*cos(alpha)))
    z4= [-r]*N
    z5=[]
    for i in range(len(s[5])):
        z5.append(-(r-sin(s[5][i])*r))
    z.extend((z0, z1, z2, z3, z4, z5))

    fy=[]
    fy0 = [element * tsk*r for element in y[0]]  # defines the function to be integrated = t_skin*y dy which becomes t_
    # skin*r*sin(theta) dtheta
    fy1 = [element * tsp for element in y[1]]
    fy2 = [element * tsk for element in y[2]]
    fy3 = [element * tsk for element in y[3]]
    fy4 = [element * tsp for element in y[4]]
    fy5 = [element * tsk*r for element in y[5]]
    fy.extend((fy0,fy1,fy2,fy3,fy4,fy5))

    fz=[]
    fz0 = [element * tsk*r for element in z[0]]  # defines the function to be integrated = t_skin*y dy which becomes
    # t_skin*r*sin(theta) dtheta
    fz1 = [element * tsp for element in z[1]]
    fz2 = [element * tsk for element in z[2]]
    fz3 = [element * tsk for element in z[3]]
    fz4 = [element * tsp for element in z[4]]
    fz5 = [element * tsk*r for element in z[5]]
    fz.extend((fz0,fz1,fz2,fz3,fz4,fz5))
    return s,z, y,fy,fz

def findboom():
    bindxlst = []
    boomlst = [[1], [], [2, 3, 4, 5], [6, 7, 8, 9], [], [10]]
    for i in range(len(s)):  # for every section
        bindxlst_i=[]
        if i == 0 or i == 2 or i == 3 or i == 5: #for sections with booms
            for k in boomlst[i]: #for every boom in the section
                err = []
                for j in range(len(s[i])):  # for all the points in the section
                    err.append((abs(boomLoc[k][1] - y[i][j])))
                bindxlst_i.append(err.index(min(err)))
        else: #for sections without booms
            bindxlst_i.append("")
        bindxlst.append(bindxlst_i)
    return bindxlst

def trapezoid(f,x):  # f should be a list with all the values of your function at a given coordinate s[i], s should be
    # a list with the coordinates in the given integration interval
    integral = 0
    intlist = []
    M = len(f)
    dx = x[1] - x[0]
    for i in range(1, M):  # The higher the N, the higher the resolution of the integral
        integral += (f[i - 1] + f[i]) * dx / 2
        intlist.append(integral)  # The i-th element in intlist represents the value of the integral at position i -->
        # the last element of the list that "trapezoid" generates is the value of the integral
    return intlist

##################### SHEAR CENTRE ######################################################################

def getsc():
    intq0r = trapezoid(qy[0], s[0])[-1] * r*r
    q2d=[]
    q3d=[]
    d = sin(alpha)*(Ca-r)
    for i in range(len(qy[2])):
        q2d.append(d* qy[2][i])
        q3d.append(d* qy[3][i])
    intq2d= trapezoid(q2d, s[2])[-1]
    intq3d= trapezoid(q3d, s[3])[-1]
    intq5d = trapezoid(qy[5], s[5])[-1] * r*r
    etha =(+intq0r + intq2d + intq3d + intq5d)
    return etha

######################## SHEAR FLOW DUE TO Y #####################################################

def Ygetqbs():
    qbs = []
    boomlst = [[1], [], [2, 3, 4, 5], [6, 7, 8, 9], [], [10]]
    for i in range(len(s)):  # for all the 6 segments
        qbs_i = []
        if i == 0 or i == 2 or i == 3 or i == 5: #for sections with booms
            m = bindxlst[i][0]
            shearjump = 0
            for j in range(len(s[i])-1):  # for all the discretized points in the segment
                if j < m:
                    qbs_i.append(-Vy / I_zz * trapezoid(fy[i], s[i])[j] + shearjump)
                elif j == m:
                    qbs_i.append(-Vy / I_zz * (boomArea * boomLoc[boomlst[i][bindxlst[i].index(m)]][1]) + qbs_i[-1])
                    shearjump += (-Vy / I_zz * (boomArea * boomLoc[boomlst[i][bindxlst[i].index(m)]][1]))
                elif j > m:
                    qbs_i.append((-Vy / I_zz * trapezoid(fy[i], s[i])[j]) + shearjump)  #  qbsi[-1] adds the constant
            # value of jump in shear due to the boom which is equal to the shear flow at the boom minus the one before
                    if bindxlst[i].index(m) != (len(bindxlst[i]) - 1): # if you didnt reach the last boom in the section yet
                        m = bindxlst[i][bindxlst[i].index(m) + 1]  # go to next boom --> increase index m
        else: #for sections without booms
            for j in range(len(s[i])-1):
                qbs_i.append(-Vy / I_zz * trapezoid(fy[i], s[i])[j])
        qbs.append(qbs_i)
    print("qbs[2] no contr=", qbs[2][-1])
    print("qbs[3] no contr=", qbs[3][-1])
    print("qbs[5] no contr=", qbs[5][-1])
    qbs[2]=[element + (qbs[0][-1]+qbs[1][-1]) for element in qbs[2]] # flow in section 1 and 2 is added to flow
    # computed for section 3
    qbs[3]=[element + qbs[2][-1] for element in qbs[3]] # flow in section 3 is added to flow computed for section 4
    qbs[5]=[element + (qbs[3][-1]-qbs[4][-1]) for element in qbs[5]] # flow in section 4 is added and flow in 5
    # subtracted to flow computed for section 6
    return qbs

def Ygetqs0():  # calculates redundant shear flow qs01, qs02 for the two cells
    intqbs = []  # calculates integral of base shear flow over each section
    for i in range(len(s)):
        intqbs.append(trapezoid(qbs[i], s[i])[-1])
        if i == 0 or i == 5:
            intqbs[-1] *= r
        if i == 4:
            intqbs[-1] *= -1
    c1 = (1 / tsk) * intqbs[0] - (1 / tsp) * intqbs[1] - 1 / tsp * intqbs[4] + 1/ tsk * intqbs[5]
    c2 = 1 / tsp * intqbs[1] + 1 / tsk * intqbs[2] + 1 / tsk * intqbs[3] + 1 / tsp * intqbs[4]
    A = array([[r * pi / tsk + 2 * r / tsp, -2 * r / tsp], [-2 * r/ tsp, 2 * ladiag / tsk + 2 * r / tsp]])
    b = array([[c1], [c2]])
    qs0 = dot(linalg.inv(A), -b)
    qs01 = qs0[0][0]
    qs02 = qs0[1][0]
    return qs01,qs02

def Ygetq():
    q=[]
    for i in range(len(s)):
        q_i=[]
        if i==0 or i==5:
            for j in range(len(qbs[i])):
                q_i.append(qbs[i][j]+qs01)
        if i==1:
            for j in range(len(qbs[i])):
                q_i.append(qbs[i][j]-qs01+qs02)
        if i==2 or i==3:
            for j in range(len(qbs[i])):
                q_i.append(qbs[i][j]+qs02)
        if i==4:
            for j in range(len(qbs[i])):
                q_i.append(qbs[i][j]+qs01-qs02)
        q.append(q_i)
    return q

######################## SHEAR FLOW DUE TO Z #####################################################

def Zgetq():
    qz = []
    boomlst = [[1], [], [2, 3, 4, 5], [6, 7, 8, 9], [], [10]]
    for i in range(len(s)):  # for all the 6 segments
        q_i = []
        if i == 0 or i == 2 or i == 3 or i == 5: #for sections with booms
            m = bindxlst[i][0]
            shearjump = 0
            for j in range(len(s[i])-1):  # for all the discretized points in the segment
                if j < m:
                    q_i.append(-Vx / I_yy* trapezoid(fz[i], s[i])[j] + shearjump)
                elif j == m:
                    q_i.append(-Vx / I_yy * (boomArea * boomLoc[boomlst[i][bindxlst[i].index(m)]][0]) + q_i[-1])
                    shearjump += (-Vx / I_yy * (boomArea * boomLoc[boomlst[i][bindxlst[i].index(m)]][0]))
                elif j > m:
                    q_i.append((-Vx / I_yy * trapezoid(fz[i], s[i])[j]) + shearjump)  #qbsi[-1] adds the constant value
                    # of jump in shear due to the boom which is equal to the shear flow at the boom minus the one before
                    if bindxlst[i].index(m) != (len(bindxlst[i]) - 1):  # if you didnt reach the last boom in the
                        # section yet
                        m = bindxlst[i][bindxlst[i].index(m) + 1]  # go to next boom --> increase index m
        else: #for sections without booms
            for j in range(len(s[i])-1):
                q_i.append(-Vx / I_yy * trapezoid(fz[i], s[i])[j])
        qz.append(q_i)
    qz[2]=[element + (qz[0][-1]+qz[1][-1]) for element in qz[2]] # flow in section 1 and 2 is added to flow computed
    # for section 3
    qz[3]=[element + qz[2][-1] for element in qz[3]] # flow in section 3 is added to flow computed for section 4
    qz[5]=[element + (qz[3][-1]-qz[4][-1]) for element in qz[5]] # flow in section 4 is added and flow in 5 subtracted
    # to flow computed for section 6
    return qz


############################# PROGRAM, SHEAR CENTRE AND SHEAR FLOWS ###################################################

# initialize values

Vy=1
Vx=1
s,z,y,fy,fz = getszyf()
bindxlst = findboom()

# determine base, redundant, final shear flows due to Vy
qbs = Ygetqbs()
qs01,qs02 = Ygetqs0()
qy= Ygetq()

# determine shear centre
zsc= getsc()-r

# determine shear flow due to Z
qz = Zgetq()

# determine shear flow due to T


##################################### DEBUG


print(trapezoid(fy[1], s[1]))


print("N=",N)
print()
print("shear centre z=", zsc)
# print()
# print("qbs[0]=",qbs[0][-1])
# print("qbs[1]=",qbs[1][-1])
# print("qbs[2]=",qbs[2][-1])
# print("qbs[3]=",qbs[3][-1])
# print("qbs[4]=",qbs[4][-1])
# print("qbs[5]=",qbs[5][-1])
# print()
# print("q[0]=",q[0][-1])
# print("q[1]=",q[1][-1])
# print("q[2]=",q[2][-1])
# print("q[3]=",q[3][-1])
# print("q[4]=",q[4][-1])
# print("q[5]=",q[5][-1])
# print("qs01=",qs01,"qs02=",qs02)


