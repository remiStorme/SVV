from numpy import *
from math import pi
from MOI import *

######################## Part I - parameters as in assignment #######################################

aircraft = "F100" # Write either A320, F100, CRJ700 or Do228 (bear in mind capitals); this is used for aerodynamic loading
Ca = 0.505  # m
la = 1.611  # m
x1 = 0.125  # m
x2 = 0.498  # m
x3 = 1.494  # m
xa = 0.25   # m
ha = 0.16  # m
tsk = 1.1/1000  # m
tsp = 2.4/1000  # m
tst = 1.2/1000  # m
hst = 13./1000   # m
wst = 17./1000   # m
nst = 11  # -
d1 = 0.00389  # m
d3 = 0.01245  # m
theta = radians(30)  # rad
P = 49.2*1000  # N
r=ha/2 #m
ladiag=sqrt((Ca-r)**2+(r)**2)

boomLoc, boomArea, z, y= centroid()
I_zz= inertiaZZ(boomArea)
I_yy= inertiaYY(boomArea, z)
# thickness=[t_skin, t_spar, t_skin, t_skin, t_spar,t_skin]

def trapezoid(f,s): # f should be a list with all the values of your function at a given coordinate s[i], s should be a list with the coordinates in the given integration interval
    integral=0
    N=len(s)
    ds=s[1]-s[0]
    for i in range (1,N): # The higher the N, the higher the resolution of the integral
        integral += (f[i-1]+f[i]) * ds / 2
    return integral



s1=linspace(0,pi/2,100000) # (start,stop,number of steps). s is the angle going from 0 to 90deg
s2=linspace(0,r,100000)
s3=linspace(0,ladiag,100000)
s4=s3
s5=s2
s6=s1
y1=[]
for i in range(len(s1)):
    y1_i=r*sin(s1[i])
    y1.append(y1_i)
y2=s2
y3=[]
for i in range(len(s3)):
    y3_i=(ladiag-s3[i])*sin(theta)
    y3.append(y3_i)
y4= [element * -1 for element in y3]
y5= [element * -1 for element in y2]
y6= [element * -1 for element in y1]


f1=[element * tsk for element in y1] # defines the function to be integrated = t_skin*y dy which becomes t_skin*r*sin(theta) dtheta
f2=[element * tsp for element in y2]
f3=[element * tsk for element in y3]
f4=[element * tsk for element in y4]
f5=[element * tsp for element in y5]
f6=[element * tsk for element in y6]
# print(trapezoid(f1,s1))
# print(trapezoid(f2,s2))
# print(trapezoid(f3,s3))
# print(trapezoid(f4,s4))
# print(trapezoid(f5,s5))
# print(trapezoid(f6,s6))
# qb1= -I_zz**(-1)*(trapezoid(f1,s1)+boomArea*boomLoc[1][1])+qb6
# qb2= -I_zz**(-1)*(trapezoid(f2,s2))
# qb3= -I_zz**(-1)*(trapezoid(f3,s3)+boomArea*boomLoc[2][1]+boomArea*boomLoc[3][1]+boomArea*boomLoc[4][1]+boomArea*boomLoc[5][1])+qb1+qb2
# qb4= -I_zz**(-1)*(trapezoid(f4,s4)+boomArea*boomLoc[6][1]+boomArea*boomLoc[7][1]+boomArea*boomLoc[8][1]+boomArea*boomLoc[9][1])+qb3
# qb5= -I_zz**(-1)*(trapezoid(f5,s5))
# qb6= -I_zz**(-1)*(trapezoid(f6,s6)+boomArea*boomLoc[10][1])+qb4-qb5



