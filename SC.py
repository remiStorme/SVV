######################## IMPORTS ########################################################

from numpy import *
from math import pi
from MOI import CrossSection
# import matplotlib.pyplot as plt

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
nst = 11  #m
d1 = 0.00389  # m
d3 = 0.01245  # m
theta = radians(30)  # rad
P = 49.2 * 1000  # N
r = ha / 2  # m
ladiag = sqrt((Ca - r) ** 2 + (r) ** 2)
cs = CrossSection()
boomLoc = cs.locBooms()
boomArea = cs.boomArea
z, y = cs.centroid()
I_zz = cs.inertiaYY()
I_yy = cs.inertiaZZ()

######################## DISCRETIZATION PARAMETERS ################################################

N = 100 # !!!!!!!! do evenly discretized cross section, not like it is now
nbooms = 11 #!!!!!!!!! include boom 0 for when Sx will be introduced

####################### IMPROVEMENTS TO IMPLEMENT ################################################

# !!!!!!!!!!!! apply simmetry coz
# !!!!!!!!!! include variable Sy=1 Newton
# !!!!!!!!!!!! fix calc of qs0 such that shear flow due to Sx can then be computed too

######################## FUNCTIONS DEFINITION #####################################################

def gety():
    s=[]
    s1 = linspace(0, (pi / 2), N)  # (start,stop,number of steps). s is the angle going from 0 to boom1, so its actually theta and the integral infact is in dtheta not ds
    s2 = linspace(0, r, N)
    s3 = linspace(0, ladiag, N)
    s4 = s3
    s5 = s2
    s6 = s1
    s.extend((s1,s2,s3,s4,s5,s6))
    y=[]
    y1 = []
    for i in range(len(s[0])):
        y1_i = r * sin(s[0][i])
        y1.append(y1_i)
    y2 = s[1]
    y3 = []
    for i in range(len(s[2])):
        y3_i = ((ladiag - s[2][i])/ladiag) *r
        y3.append(y3_i)
    y4 = [element-r for element in y3]
    y5 = [element * -1 for element in y2]
    y6=[]
    for i in range(len(y1)):
        y6.append(-y1[N-1-i])
    y.extend((y1,y2,y3,y4,y5,y6))
    return s, y

def findboom(ysec, nboom): #to find the position of boom 1, section=s[0], nboom=1
    err = []
    ind = []
    for i in range(len(ysec)):
        err.append((abs(boomLoc[nboom][1] - ysec[i])))
    ind.append(err.index(min(err)))
    return ind

def boomindxlst():
    i1=findboom(y[0],1)
    i2= findboom(y[2], 2)
    i3=findboom(y[2], 3)
    i4=findboom(y[2], 4)
    i5=findboom(y[2], 5)
    i6=findboom(y[3], 6)
    i7=findboom(y[3], 7)
    i8=findboom(y[3], 8)
    i9=findboom(y[3], 9)
    i10=findboom(y[5], 10)
    bindxlst=[i1]+[[]]+[i2+i3+i4+i5]+[i6+i7+i8+i9]+[[]]+[i10]
    return bindxlst

def trapezoid(f,s):  # f should be a list with all the values of your function at a given coordinate s[i], s should be a list with the coordinates in the given integration interval
    integral = 0
    intlist = []
    M = len(s)
    ds = s[1] - s[0]
    for i in range(1, M):  # The higher the N, the higher the resolution of the integral
        integral += (f[i - 1] + f[i]) * ds / 2
        intlist.append(integral) # The i-th element in intlist represents the value of the integral at position i --> the last element of the list that "trapezoid" generates is the value of the integral
    return intlist

def getf():
    f=[]
    f1 = [element * tsk*r for element in y[0]]  # defines the function to be integrated = t_skin*y dy which becomes t_skin*r*sin(theta) dtheta
    f2 = [element * tsp*r for element in y[1]]
    f3 = [element * tsk*r for element in y[2]]
    f4 = [element * tsk*r for element in y[3]]
    f5 = [element * tsp*r for element in y[4]]
    f6 = [element * tsk*r for element in y[5]]
    f.extend((f1,f2,f3,f4,f5,f6))
    return f

def getqbs():
    qbs=[]
    boomlst = [[1], [], [2, 3, 4, 5], [6, 7, 8, 9], [], [10]]
    for i in range(len(s)): #for all the 6 segments
        qbs_i=[]
        if i==0 or i==2 or i==3 or i==5:
            m=bindxlst[i][0]
            shearjump=0
            for j in range(len(s[i])): #for all the discretized points in the segment
                if j<m:
                    # print("i=",i,"j=",j,"m=",m)
                    qbs_i.append(-1/I_zz * trapezoid(f[i],s[i])[j]+shearjump)
                elif j==m:
                    # print("i=", i, "j=", j, "m=", m)
                    qbs_i.append(-1/I_zz * (boomArea*boomLoc[boomlst[i][bindxlst[i].index(m)]][1])+qbs_i[-1])
                    # print("started from the bottom now we here: boom number", boomlst[i][bindxlst[i].index(m)],"m=", m)
                    shearjump += (-1 / I_zz * (boomArea * boomLoc[boomlst[i][bindxlst[i].index(m)]][1]))
                    # print("shearjump", shearjump)
                elif j>m:
                    # print("i=", i, "j=", j, "m=", m)
                    qbs_i.append((-1/I_zz * trapezoid(f[i],s[i])[j-1])+shearjump) # j-1 because trapezoid() leads to a list of values of areas that sum up to an integral, so number of elements in trapezoid is always 1 less than integrated function. qbsi[-1] adds the constant value of jump in shear due to the boom which is equal to the shear flow at the boom minus the one before
                    if bindxlst[i].index(m) != (len(bindxlst[i]) - 1): # if you didnt reach the last boom in the section yet
                        m = bindxlst[i][bindxlst[i].index(m) + 1] # go to next boom --> increase index m
                # print("SHEAR FLOW",qbs_i[-1])
        else:
            for j in range(len(s[i])):
                # print("i=", i, "j=", j, "m=", m)
                qbs_i.append(-1 * I_zz * trapezoid(f[i], s[i])[j-1])
                # print("SHEAR FLOW", qbs_i[-1])
        qbs.append(qbs_i)
    qbs[2]=[element + (qbs[0][-1]+qbs[1][-1]) for element in qbs[2]] # flow in section 1 and 2 is added to flow computed for section 3
    qbs[3]=[element + qbs[2][-1] for element in qbs[3]] # flow in section 3 is added to flow computed for section 4
    qbs[5]=[element + (qbs[3][-1]-qbs[4][-1]) for element in qbs[5]] # flow in section 4 is added and flow in 5 subtracted to flow computed for section 6
    return qbs

def getintqbs():
    intqbs=[]
    for k in range(len(s)):
        intqbs.append(trapezoid(qbs[k],s[k])[-1])
    return intqbs
############################# PROGRAM ###################################################

s,y=gety()
bindxlst=boomindxlst()
f=getf()
qbs=getqbs()
intqbs=getintqbs()
print(intqbs)
c1=-(pi/(2*tsk))*intqbs[0]-(1/tsp)*intqbs[1]-1/tsp*intqbs[4]-pi/(2*tsk)*intqbs[5]

c2=-r/tsp*intqbs[1]-ladiag/tsk*intqbs[2]-ladiag/tsk*intqbs[3]-r/tsp*intqbs[4]
print(c1,c2)


########################### END PROGRAM #################################################















# def getqs0():
#     num_i=[]
#     den_i=[]
#     for i in range(1,5):
#         num_i.append(trapezoid(qbs[i],s[i])[-1]) # the last element of the list that "trapezoid" generates is the value of the integral
#         den_i.append(trapezoid(N*[1],s[i])[-1])
#     num=sum(num_i)
#     den=sum(den_i)
#     qs0 = -num / den
#     return qs0,den