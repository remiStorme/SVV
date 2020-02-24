from numpy import *
from math import pi
from MOI import CrossSection

######################## GEOMETRY F100 #######################################


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

#############################################################################
discpts = 1000
nbooms = 11

def gety():
    s=[]
    s1 = linspace(0, (pi / 2), discpts)  # (start,stop,number of steps). s is the angle going from 0 to boom1
    s2 = linspace(0, r, discpts)
    s3 = linspace(0, ladiag, discpts)
    s4 = s3
    s5 = s2
    s6 = s1
    s.extend((s1,s2,s3,s4,s5,s6))
    y=[]
    y1 = []
    for i in range(len(s1)):
        y1_i = r * sin(s1[i])
        y1.append(y1_i)
    y2 = s2
    y3 = []
    for i in range(len(s3)):
        y3_i = (ladiag - s3[i]) * sin(theta)
        y3.append(y3_i)
    y4 = [element * -1 for element in y3]
    y5 = [element * -1 for element in y2]
    y6 = [element * -1 for element in y1]
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
    N = len(s)
    ds = s[1] - s[0]
    for i in range(1, N):  # The higher the N, the higher the resolution of the integral
        integral += (f[i - 1] + f[i]) * ds / 2
        intlist.append(integral)
    return intlist

def getf():
    f=[]
    f1 = [element * tsk for element in y[0]]  # defines the function to be integrated = t_skin*y dy which becomes t_skin*r*sin(theta) dtheta
    f2 = [element * tsp for element in y[1]]
    f3 = [element * tsk for element in y[2]]
    f4 = [element * tsk for element in y[3]]
    f5 = [element * tsp for element in y[4]]
    f6 = [element * tsk for element in y[5]]
    f.extend((f1,f2,f3,f4,f5,f6))
    return f

def getqbs():
    qbs=[]
    qbs_i=[]
    boomlst = [[1], [], [2, 3, 4, 5], [6, 7, 8, 9], [], [10]]
    for i in range(len(s)): #for all the 6 segments
        for j in range(len(s[i])): #for all the discretized points in the segment
           for m in (bindxlst[i]): #for all the booms in the segment
                if j<m:
                    print("i=",i,"j=",j,"m=",m) #check
                    qbs_i.append(-1 * I_zz * trapezoid(f[i],s[i])[j])
                elif j==m:
                    qbs_i.append(-1 * I_zz * boomArea*boomLoc[boomlst[i][bindxlst[i].index(m)]][1])
                    print("started from the bottom now we here") #check
                elif j>m:
                    qbs_i.append(-1 * I_zz * trapezoid(f[i],s[i])[j-1]) #why j-1 find out
                    print("i=",i,"j=",j,"m=",m)
                else:
                    print("uhm wuuuuut lorenza get it togetha")
        qbs.append(qbs_i)
    return qbs

############################# PROGRAM ###################################################

s,y=gety()
bindxlst=boomindxlst()
print(bindxlst)
f=getf()
qbs=getqbs()
print(len(qbs))


########################### END PROGRAM #################################################



#     qb=[]
#     qb1 = -I_zz ** (-1) * (int[0] + boomArea * boomLoc[1][1])
#     qb2 = -I_zz ** (-1) * (int[1])
#     qb3 = -I_zz ** (-1) * (int[2] + boomArea * boomLoc[2][1] + boomArea * boomLoc[3][1] + boomArea * boomLoc[4][1] + boomArea *
#                 boomLoc[5][1]) + qb1 + qb2
#     qb4 = -I_zz ** (-1) * (int[3] + boomArea * boomLoc[6][1] + boomArea * boomLoc[7][1] + boomArea * boomLoc[8][1] + boomArea *
#                 boomLoc[9][1]) + qb3
#     qb5 = -I_zz ** (-1) * (int[4])
#     qb6 = -I_zz ** (-1) * (int[5] + boomArea * boomLoc[10][1]) + qb4 - qb5
#     qb.extend((qb1,qb2,qb3,qb4,qb5,qb6))
#     return qb

