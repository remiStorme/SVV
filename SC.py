######################## IMPORTS ########################################################
import itertools
from numpy import *
from math import pi
from MOI import CrossSection
import Internal_Forces_and_Deflections as ifd
from matplotlib import pyplot as plt
######################## GEOMETRY F100 ########################################################

Ca = 0.505  # m
la = 1.611  # m
x1 = 0.125  # m
x2 = 0.498  # m
x3 = 1.494  # m
xa = 0.25  # m
ha = 0.161  # m
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
l = sqrt(r ** 2 + (Ca - r) ** 2)
ladiag = sqrt((Ca - r) ** 2 + (r) ** 2)
alpha = arctan(r/(Ca-r))
A_i = pi * r ** 2 / 2
A_ii = ha * (Ca - r) / 2
Bi = wst * tst + hst * tst
cs = CrossSection()
boomLoc = cs.locBooms()
boomArea = 0#3.6 * 10 ** (-5)
I_zz = cs.inertiaZZ()
I_yy = cs.inertiaYY()
zcentr,ycentr=0.20362591085157106,0.0 #taken from verification model because numerical model Zsc is wrong
G = 28e9

######################## DISCRETIZATION PARAMETERS ################################################

N = 100
nbooms = 11

######################## INTEGRATION TOOLS #####################################################

def getszyf():
    s = []
    s0 = linspace(0, (pi / 2),N)  # (start,stop,number of steps). s is the angle going from 0 to boom1, so its actually theta and the integral infact is in dtheta not ds
    s1 = linspace(0, r, N)
    s2 = linspace(0, ladiag, N)
    s3 = s2
    s4 = linspace(0, -r, N)
    s5 = linspace((-pi/2), 0,N)
    s.extend((s0, s1, s2, s3, s4, s5))

    y = []
    y0 = []
    for i in range(len(s[0])):
        y0.append(r * sin(s[0][i]))
    y1 = s[1]
    y2 = []
    for i in range(len(s[2])):
        y2_i = r - (s[2][i])*sin(alpha)
        y2.append(y2_i)
    y3 = []
    for i in range(len(s[2])):
        y3_i = - (s[2][i])*sin(alpha)
        y3.append(y3_i)
    y4 = [element * -1 for element in y1]
    y5=  []
    for i in range(len(s[5])):
        y5.append(-(cos(s[0][i])*r))
    y.extend((y0, y1, y2, y3, y4, y5))

    z= []
    z0=[]
    for i in range(len(s[0])):
        z0.append(((cos(s[0][i])*r)+(zcentr-r)))
    z1= [zcentr-r]*N
    z2=[]
    for i in range(len(s[2])):
        z2_i = (zcentr-r)-s[2][i]*cos(alpha)
        z2.append(z2_i)
    z3=[]
    for i in range(len(s[3])):
        z3_i = -(Ca-zcentr)+s[2][i]*cos(alpha)
        z3.append(z3_i)
    z4= [zcentr-r]*N
    z5=[]
    for i in range(len(s[5])):
        z5.append(-(-(sin(s[0][i])*r+(zcentr-r))))
    z.extend((z0, z1, z2, z3, z4, z5))
    fy=[]
    fy0 = [element * tsk*r for element in y[0]]  # defines the function to be integrated = t_skin*y dy which becomes t_skin*r*sin(theta) dtheta
    fy1 = [element * tsp for element in y[1]]
    fy2 = [element * tsk for element in y[2]]
    fy3 = [element * tsk for element in y[3]]
    fy4 = [element * tsp for element in y[4]]
    fy5 = [element * tsk*r for element in y[5]]
    fy.extend((fy0,fy1,fy2,fy3,fy4,fy5))

    fz=[]
    fz0 = [element * tsk*r for element in z[0]]  # defines the function to be integrated = t_skin*y dy which becomes t_skin*r*sin(theta) dtheta
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

def trapezoid(f,x):  # f should be a list with all the values of your function at a given coordinate s[i], s should be a list with the coordinates in the given integration interval
    integral = 0
    intlist = []
    M = len(f)
    dx = x[1] - x[0]
    for i in range(1, M):  # The higher the N, the higher the resolution of the integral
        integral += (f[i - 1] + f[i]) * dx / 2
        intlist.append(integral)  # The i-th element in intlist represents the value of the integral at position i --> the last element of the list that "trapezoid" generates is the value of the integral
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

def Ygetqbs(Vy):
    qbs = []
    boomlst = [[1], [], [2, 3, 4, 5], [6, 7, 8, 9], [], [10]]
    for i in range(len(s)):  # for all the 6 segments
        qbs_i = []
        if i == 0 or i == 2 or i==3 or i==5: #for sections with booms
            m = bindxlst[i][0]
            shearjump = 0
            for j in range(len(s[i])-1):  # for all the discretized points in the segment
                if j < m:
                    qbs_i.append(-Vy / I_zz * trapezoid(fy[i], s[i])[j] + shearjump)
                elif j == m:
                    qbs_i.append(-Vy / I_zz * (boomArea * boomLoc[boomlst[i][bindxlst[i].index(m)]][1]) + qbs_i[-1])
                    shearjump += (-Vy / I_zz * (boomArea * boomLoc[boomlst[i][bindxlst[i].index(m)]][1]))
                elif j > m:
                    qbs_i.append((-Vy / I_zz * trapezoid(fy[i], s[i])[j]) + shearjump)  #  qbsi[-1] adds the constant value of jump in shear due to the boom which is equal to the shear flow at the boom minus the one before
                    if bindxlst[i].index(m) != (len(bindxlst[i]) - 1):  # if you didnt reach the last boom in the section yet
                        m = bindxlst[i][bindxlst[i].index(m) + 1]  # go to next boom --> increase index m
        else:
            for j in range(len(s[i]) - 1):
                qbs_i.append(-Vy / I_zz * trapezoid(fy[i], s[i])[j])
        qbs.append(qbs_i)
    qbs[2]=[element + ((qbs[0][-1])+(qbs[1][-1])) for element in qbs[2]] # flow in section 1 and 2 is added to flow computed for section 3
    qbs[3]=[element + (qbs[2][-1]) for element in qbs[3]] # flow in section 3 is added to flow computed for section 4
    qbs[5]=[element + ((qbs[3][-1])-(qbs[4][-1])) for element in qbs[5]] # flow in section 4 is added and flow in 5 subtracted to flow computed for section 6
    return qbs

def Ygetqs0(qbs):  # calculates redundant shear flow qs01, qs02 for the two cells
    intqbs = []  # calculates integral of base shear flow over each section
    for i in range(len(s)):
        intqbs.append(trapezoid(qbs[i], s[i])[-1])
        if i == 0 or i == 5:
            intqbs[-1] *= r
    c1 = (1 / tsk) * intqbs[0] - (1 / tsp) * intqbs[1] - 1 / tsp * intqbs[4] + 1/ tsk * intqbs[5]
    c2 = 1 / tsp * intqbs[1] + 1 / tsk * intqbs[2] + 1 / tsk * intqbs[3] + 1 / tsp * intqbs[4]
    A = array([[r * pi / tsk + 2 * r / tsp, -2 * r / tsp], [-2 * r/ tsp, 2 * ladiag / tsk + 2 * r / tsp]])
    b = array([[c1], [c2]])
    qs0 = dot(linalg.inv(A), -b)
    qs01 = qs0[0][0]
    qs02 = qs0[1][0]
    return qs01,qs02

def Ygetq(qbs,qs01,qs02):
    q=[]
    for i in range(len(s)):
        q_i=[]
        if i==0:
            for j in range(len(qbs[i])):
                q_i.append(qbs[i][j]+qs01)
        if i==1:
            for j in range(len(qbs[i])):
                q_i.append(qbs[i][j]-qs01+qs02)
        if i==2:
            for j in range(len(qbs[i])):
                q_i.append(qbs[i][j]+qs02)
        if i==3:
            for j in range(len(qbs[i])):
                q_i.append(qbs[i][j]+qs02)
        if i==4:
            for j in range(len(qbs[i])):
                q_i.append(qbs[i][j]-qs01+qs02)
        if i==5:
            for j in range(len(qbs[i])):
                q_i.append(qbs[i][j]+qs01)
        q.append(q_i)
    return q

######################## SHEAR FLOW DUE TO Z #####################################################

def Zgetq(Vz):
    qz = []
    boomlst = [[1], [], [2, 3, 4, 5], [6, 7, 8, 9], [], [10]]
    for i in range(len(s)):  # for all the 6 segments
        q_i = []
        if i == 0 or i == 2: #for sections with booms
            m = bindxlst[i][0]
            shearjump = 0
            for j in range(len(s[i])-1):  # for all the discretized points in the segment
                if j < m:
                    q_i.append(-Vz / I_yy * trapezoid(fz[i], s[i])[j] + shearjump)
                elif j == m:
                    q_i.append(-Vz / I_yy * (boomArea * -(boomLoc[boomlst[i][bindxlst[i].index(m)]][0]-zcentr)) + q_i[-1])
                    shearjump += (-Vz / I_yy * (boomArea * (-(boomLoc[boomlst[i][bindxlst[i].index(m)]][0]-zcentr))))
                elif j > m:
                    q_i.append((-Vz / I_yy * trapezoid(fz[i], s[i])[j]) + shearjump)  #  qbsi[-1] adds the constant value of jump in shear due to the boom which is equal to the shear flow at the boom minus the one before
                    if bindxlst[i].index(m) != (len(bindxlst[i]) - 1):  # if you didnt reach the last boom in the section yet
                        m = bindxlst[i][bindxlst[i].index(m) + 1]  # go to next boom --> increase index m
        if i==1: #for spar without booms
            for j in range(len(s[i])-1):
                q_i.append(-Vz / I_yy * trapezoid(fz[i], s[i])[j])
        if i==3:
            qz[2] = [element + ((qz[0][-1]) + (qz[1][-1])) for element in qz[2]]  # flow in section 1 and 2 is added to flow computed for section 3
            q_i=[-element for element in qz[2][::-1]]
        if i==4:
            q_i=[-element for element in qz[1]]
        if i==5:
            q_i=[-element for element in qz[0][::-1]]
        qz.append(q_i)
    return qz

######################## SHEAR FLOW DUE TO T AND TORSIONAL CONSTANT #####################################################

def Tgetq(T):
    A = zeros((3,3))
    b = zeros(3)

    ### First row
    A[0, 0] = 2. * A_i
    A[0, 1] = 2. * A_ii
    b[0] = -T

    ### Second row
    A[1, 0] = (r * pi / tsk + 2 * r / tsp) / (2 * A_i)
    A[1, 1] = (-2 * r / tsp) / (2 * A_i)
    A[1, 2] = -1.
    b[1] = 0.

    ### Third row
    A[2, 0] = (-2 * r / tsp) / (2 * A_ii)
    A[2, 1] = (2 * ladiag / tsk + 2 * r / tsp) / (2 * A_ii)
    A[2, 2] = -1
    b[2] = 0.

    F = linalg.solve(A, b)
    q01n = F[0]
    q02n = F[1]
    J = -1/F[2]

    Tq=[]
    for i in range(len(s)):
        Tq_j = []
        if i==0 or i==5:
            for j in range(len(s[i])-1):
                Tq_j.append(q01n)
        if i==1:
            for j in range(len(s[i])-1):
                Tq_j.append(-q01n+q02n)
        if i==2 or i==3:
            for j in range(len(s[i])-1):
                Tq_j.append(q02n)
        if i==4:
            for j in range(len(s[i])-1):
                Tq_j.append(-q01n+q02n)
        Tq.append(Tq_j)
    return Tq, J

############################# PROGRAM ###################################################

# initialize values
s,z,y,fy,fz = getszyf()
bindxlst = findboom()
x_T,  y_T =ifd.x_T,  ifd.y_T #x_T is the slice at which Tx is being evaluated and y_sy istotal torque on x axis around shear centre generated by ae load and reaction forces
x_sy, y_sy = ifd.x_sy, ifd.y_sy #x_sy is the slice at which sy is being evaluated and y_sy is the magnitude of Sy at that slice
x_sz, y_sz = ifd.x_sz, ifd.y_sz #x_sz is the slice at which sz is being evaluated and y_sz is the magnitude of Sz at that slice

# determine base, redundant, final shear flows due to Vy
qy=[]
for i in range(len(x_sy)):
    qy_slice_i=[]
    # Vy=y_sy[i][0]
    Vy=1                                                 # change here!!
    qbs_slice_i=Ygetqbs(Vy)
    qs01_slice_i,qs02_slice_i = Ygetqs0(qbs_slice_i)
    qy_slice_i=(Ygetq(qbs_slice_i,qs01_slice_i,qs02_slice_i))
    qy.append(qy_slice_i)

# determine shear flow due to Vz
qz=[]
for i in range(len(x_sz)): #for 40 slices
    qz_slice_i=[]
    # Vz=y_sz[i][0]
    Vz=0                                                # change here!!
    qz_slice_i=(Zgetq(Vz))
    qz.append(qz_slice_i)

# determine shear flow due to T
qt=[]
for i in range(len(x_T)):
    qt_slice_i=[]
    # T=y_T[i][0]
    T=1                                                 # change here!!
    qt_slice_i,J=(Tgetq(T))
    qt.append(qt_slice_i)

# determine maximum shear stress along aileron

shearstress=[]
for i in range(len(qy)): #for all 40 slices
    qcombo_j=[]
    for j in range(len(qy[0])): #for all sections in one slice
        qcombo_k=[]
        for k in range(len(qy[0][0])-1): #for all points in one section in one slice
            if j==0 or j==2 or j==3 or j==5:
                qcombo_k.append((qy[i][j][k]+qz[i][j][k]+qt[i][j][k])/tsk)
                maxshstrsec = max(qcombo_k)
            if j == 1 or j == 4:
                qcombo_k.append((qy[i][j][k] + qz[i][j][k] + qt[i][j][k])/tsp)
                maxshstrsec= max(qcombo_k)
        qcombo_j.append(maxshstrsec)
        maxshstrslice=max(qcombo_j)
    shearstress.append(maxshstrslice)
maxss=max(shearstress)
Nslicemaxss=shearstress.index(max(shearstress))
maxsslocation=(la/40)*Nslicemaxss
print("Maximum shear stress=",maxss/10**6,"[MPa]","at location x=",maxsslocation,"[m] (from the fuselage)")

# determine total shear flow for crossection # Ncs=0
Ncs=0
qy = array((qy[Ncs])).flatten()
qz = array((qz[Ncs])).flatten()
qt = array((qt[Ncs])).flatten()

q = qy + qz + qt

# determine shear stress in the crossection # Ncs=0
zvaluesregion1 = (array(z[0])[:-1] + array(z[0])[1:]) / 2
zvaluesregion2 = (array(z[1])[:-1] + array(z[1])[1:]) / 2
zvaluesregion3 = (array(z[2])[:-1] + array(z[2])[1:]) / 2
zvaluesregion4 = (array(z[3])[:-1] + array(z[3])[1:]) / 2
zvaluesregion5 = (array(z[4])[:-1] + array(z[4])[1:]) / 2
zvaluesregion6 = (array(z[5])[:-1] + array(z[5])[1:]) / 2

yvaluesregion1 = (array(y[0])[:-1] + array(y[0])[1:]) / 2
yvaluesregion2 = (array(y[1])[:-1] + array(y[1])[1:]) / 2
yvaluesregion3 = (array(y[2])[:-1] + array(y[2])[1:]) / 2
yvaluesregion4 = (array(y[3])[:-1] + array(y[3])[1:]) / 2
yvaluesregion5 = (array(y[4])[:-1] + array(y[4])[1:]) / 2
yvaluesregion6 = (array(y[5])[:-1] + array(y[5])[1:]) / 2

tau0 = q[:N-1] / tsk
tau1 = q[N-1:2*(N-1)] / tsp
tau2 = q[2*(N-1):3*(N-1)] / tsk
tau3 = q[3*(N-1):4*(N-1)] / tsk
tau4 = q[4*(N-1):5*(N-1)] / tsp
tau5 = q[5*(N-1):] / tsk

z = concatenate((zvaluesregion1,zvaluesregion2,zvaluesregion3,zvaluesregion4,zvaluesregion5,zvaluesregion6))
y = concatenate((yvaluesregion1,yvaluesregion2,yvaluesregion3,yvaluesregion4,yvaluesregion5,yvaluesregion6))
tau = concatenate((tau0,tau1,tau2,tau3,tau4,tau5))

# Von Mises and Bending Stress sigma from selected My, Mz, tau for crossection # Ncs=0
My = 0                                                            # change here!!                                                          # change here!!
sigma = My * z / I_yy + T * y / I_zz
VM = sqrt(sigma**2 + 3 * tau**2)

# determine torsional stiffness
Tq,J= Tgetq(1)
print("Torsional stiffness=", J,"[Nm/rad]")

# determine shear centre
qbs=Ygetqbs(1)
qs01,qs02=Ygetqs0(qbs)
qy=Ygetq(qbs,qs01,qs02)
zsc= getsc()-r
print("Shear centre location=",zsc,"[m] from the spar, right is positive")

# plot shear flow
fig, ax = plt.subplots()
sfplot=plt.scatter(z,y, c=q,alpha=1,cmap='jet',vmin = -max(abs(q)),vmax = max(abs(q)))
ax.set_xlim(0.3,-0.35,)
cbar = plt.colorbar(sfplot)
cbar.set_label("Shear flow [N/m]")
plt.xlabel("z [m]")
plt.ylabel("y [m]")
plt.show()

# plot shear stress
fig, ax = plt.subplots()
sfplot=plt.scatter(z,y, c=tau,alpha=1,cmap='jet',vmin = -max(abs(tau)),vmax = max(abs(tau)))
ax.set_xlim(0.3,-0.35,)
cbar = plt.colorbar(sfplot)
cbar.set_label("Shear stress [Pa]")
plt.xlabel("z [m]")
plt.ylabel("y [m]")
plt.show()

# plot bending stress
fig, ax = plt.subplots()
sfplot=plt.scatter(z,y, c=sigma,alpha=1,cmap='jet',vmin = -max(abs(sigma)),vmax = max(abs(sigma)))
ax.set_xlim(0.3,-0.35,)
cbar = plt.colorbar(sfplot)
cbar.set_label("Bending stress [Pa]")
plt.xlabel("z [m]")
plt.ylabel("y [m]")
plt.show()

# plot Von Mises stress
fig, ax = plt.subplots()
sfplot=plt.scatter(z,y, c=VM,alpha=1,cmap='jet',vmin = 0,vmax = max(abs(VM)))
ax.set_xlim(0.3,-0.35,)
cbar = plt.colorbar(sfplot)
cbar.set_label("Von Mises [N/m]")
plt.xlabel("z [m]")
plt.ylabel("y [m]")
plt.show()

#################################################################################################### END