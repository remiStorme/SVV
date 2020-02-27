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




x_sy, y_sy = ifd.x_sy, ifd.y_sy
x_sz, y_sz = ifd.x_sz, ifd.y_szK

print("x_sy",x_sy)


torque_ae = ae.t_w #torque per cross section down the span

lift_ae = ae.w_x #resultant lift per cross section

#vector x = [Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Ar, C1, C2, C3, C4, C5]
uks = rr.uks

print(uks)

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
cs = CrossSection()
boomLoc = cs.locBooms()
boomArea = 3.456 * 10 ** (-5)
z, y = cs.centroid()
I_zz = 4.686*10**-6 #from verification
I_yy = cs.inertiaZZ()

######################## DISCRETIZATION PARAMETERS ################################################

N = 100 # !!!!!!!! do evenly discretized cross section, not like it is now
nbooms = 11 #!!!!!!!!! include boom 0 for when Sx will be introduced
Vy=1

######################## FUNCTIONS DEFINITION #####################################################

def getsyf():
    s = []
    s1 = linspace(0, (pi / 2),N)  # (start,stop,number of steps). s is the angle going from 0 to boom1, so its actually theta and the integral infact is in dtheta not ds
    s2 = linspace(0, r, N)
    s3 = linspace(0, ladiag, N)
    s4 = s3
    s5 = linspace(0, -r, N)
    s6 = linspace(-(pi / 2),0,N)
    s.extend((s1, s2, s3, s4, s5, s6))
    y = []
    y1 = []
    for i in range(len(s[0])):
        y1_i = r * sin(s[0][i])
        y1.append(y1_i)
    y2 = s[1]
    y3 = []
    for i in range(len(s[2])):
        y3_i = ((ladiag - s[2][i]) / ladiag) * r
        y3.append(y3_i)
    y4 = [element - r for element in y3]
    y5 = [element * -1 for element in y2]
    y6=  []
    for i in range(len(y1)):
        y6.append(-y1[N - 1 - i])
    y.extend((y1, y2, y3, y4, y5, y6))
    f=[]
    f1 = [element * tsk*r for element in y[0]]  # defines the function to be integrated = t_skin*y dy which becomes t_skin*r*sin(theta) dtheta
    f2 = [element * tsp for element in y[1]]
    f3 = [element * tsk for element in y[2]]
    f4 = [element * tsk for element in y[3]]
    f5 = [element * tsp for element in y[4]]
    f6 = [element * tsk*r for element in y[5]]
    f.extend((f1,f2,f3,f4,f5,f6))
    return s, y,f

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

def trapezoid(f,x):  # f should be a list with all the values of your function at a given coordinate s[i], s should be a list with the coordinates in the given integration interval
    integral = 0
    intlist = []
    M = len(f)
    dx = x[1] - x[0]
    for i in range(1, M):  # The higher the N, the higher the resolution of the integral
        integral += (f[i - 1] + f[i]) * dx / 2
        intlist.append(integral)  # The i-th element in intlist represents the value of the integral at position i --> the last element of the list that "trapezoid" generates is the value of the integral
    return intlist

def getqbs():
    qbs = []
    boomlst = [[1], [], [2, 3, 4, 5], [6, 7, 8, 9], [], [10]]
    for i in range(len(s)):  # for all the 6 segments
        qbs_i = []
        if i == 0 or i == 2 or i == 3 or i == 5: #for sections with booms
            m = bindxlst[i][0]
            shearjump = 0
            for j in range(len(s[i])-1):  # for all the discretized points in the segment
                if j < m:
                    qbs_i.append(-Vy / I_zz * trapezoid(f[i], s[i])[j] + shearjump)
                elif j == m:
                    qbs_i.append(-Vy / I_zz * (boomArea * boomLoc[boomlst[i][bindxlst[i].index(m)]][1]) + qbs_i[-1])
                    shearjump += (-Vy / I_zz * (boomArea * boomLoc[boomlst[i][bindxlst[i].index(m)]][1]))
                elif j > m:
                    qbs_i.append((-Vy / I_zz * trapezoid(f[i], s[i])[j]) + shearjump)  # j-1 because trapezoid() leads to a list of values of areas that sum up to an integral, so number of elements in trapezoid is always 1 less than integrated function. qbsi[-1] adds the constant value of jump in shear due to the boom which is equal to the shear flow at the boom minus the one before
                    if bindxlst[i].index(m) != (len(bindxlst[i]) - 1):  # if you didnt reach the last boom in the section yet
                        m = bindxlst[i][bindxlst[i].index(m) + 1]  # go to next boom --> increase index m
        else: #for sections without booms
            for j in range(len(s[i])-1):
                qbs_i.append(-Vy / I_zz * trapezoid(f[i], s[i])[j])
        qbs.append(qbs_i)
    # print(qbs[0])
    # print(qbs[1][-1])
    # print(qbs[2][-1])
    # print(qbs[3][-1])
    # print(qbs[4][-1])
    # print(qbs[5][-1])
    qbs[2]=[element + (qbs[0][-1]+qbs[1][-1]) for element in qbs[2]] # flow in section 1 and 2 is added to flow computed for section 3
    qbs[3]=[element + qbs[2][-1] for element in qbs[3]] # flow in section 3 is added to flow computed for section 4
    qbs[5]=[element + (qbs[3][-1]-qbs[4][-1]) for element in qbs[5]] # flow in section 4 is added and flow in 5 subtracted to flow computed for section 6
    # print("adding contributions of previous ")
    # print(qbs[0][-1])
    # print(qbs[1][-1])
    # print(qbs[2][-1])
    # print(qbs[3][-1])
    # print(qbs[4][-1])
    # print(qbs[5][-1])
    return qbs

def getqs0():  # calculates redundant shear flow qs01, qs02 for the two cells
    intqbs = []  # calculates integral of base shear flow over each section
    for i in range(len(s)):
        intqbs.append(trapezoid(qbs[i], s[i])[-1])
        if i == 0 or i == 5:
            intqbs[-1] *= r
        if i == 4:
            intqbs[-1] *= -1
    # c1 = (pi / (2 * tsk)) * intqbs[0] - (r / tsp) * intqbs[1] + r / tsp * intqbs[4] + pi / (2 * tsk) * intqbs[5]
    # c2 = r / tsp * intqbs[1] + ladiag / tsk * intqbs[2] + ladiag / tsk * intqbs[3] - r / tsp * intqbs[4]
    # print(intqbs)
    c1 = (1 / tsk) * intqbs[0] - (1 / tsp) * intqbs[1] - 1 / tsp * intqbs[4] + 1/ tsk * intqbs[5]
    c2 = 1 / tsp * intqbs[1] + 1 / tsk * intqbs[2] + 1 / tsk * intqbs[3] + 1 / tsp * intqbs[4]
    A = array([[r * pi / tsk + 2 * r / tsp, -2 * r / tsp], [-2 * r/ tsp, 2 * ladiag / tsk + 2 * r / tsp]])
    b = array([[c1], [c2]])
    # print(A,b)
    qs0 = dot(linalg.inv(A), -b)
    qs01 = qs0[0][0]
    qs02 = qs0[1][0]
    return qs01,qs02

def getqfinal():
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

def getsc():
    intq0r = trapezoid(q[0], s[0])[-1] * r*r
    q2d=[]
    q3d=[]
    alpha = arctan(r/(Ca-r))
    d = sin(alpha)*(Ca-r)
    for i in range(len(q[2])):
        q2d.append(d* q[2][i])
        q3d.append(d* q[3][i])
    intq2d= trapezoid(q2d, s[2])[-1]
    intq3d= trapezoid(q3d, s[3])[-1]
    intq5d = trapezoid(q[5], s[5])[-1] * r*r
    etha =(+intq0r + intq2d + intq3d + intq5d)
    return etha

############################# PROGRAM ###################################################
s,y,f = getsyf()
bindxlst = findboom()
qbs = getqbs()
qs01,qs02 = getqs0()
q=getqfinal()
zsc=getsc()-0.08
print("shear centre is at z:",zsc)
########################### END PROGRAM #################################################
