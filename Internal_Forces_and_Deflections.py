import numpy as np
# import matplotlib.pyplot as plt
import math as m
import Reaction_Forces_Aero_Load as rf
import Interpolator_Integrate_Cubic as ii
import AEload as ae

#RELEVANT VARIABLES
Ca = 0.505  # m
la = 1.611  # m
la_lim = 1.40 #m
x1 = 0.125  # m
x2 = 0.498  # m
x3 = 1.494  # m
xa = 0.245   # m
ha = 0.161  # m
r = ha/2
tsk = 1.1/1000  # m
tsp = 2.4/1000  # m
tst = 1.2/1000  # m
hst = 13./1000   # m
wst = 17./1000   # m
nst = 11  # -
d1 = 0.00389  # m
d3 = 0.01245  # m
theta = m.radians(30)  # rad
P = 49.2*1000  # N

A_t = 0.002686478946739162
yc = -4.339264140786877e-19
zc = -0.24048766835061938


#MATERIAL PROPERTIES
E       = 72.9*10**9       # E-modulus (Pa)
G       = 27.1*10**9       # G-modulus (Pa)

#HINGE AND ACTUATOR LOCATIONS: TEMPORARY VARIABLES
xaj     = x2 - xa/2
xp      = x2 + xa/2

#SECTION PROPERTIES
z_sc    = -0.08499063497059493
J       = 7.649955726444055*10**-6
I_yy    = 4.5925790464352304*10**-5
I_zz    = 4.686331664359035*10**-6


#retrieving values for 12 unknowns from the matrix system
uks = rf.uks
R_1y, R_2y, R_3y, R_1z, R_2z, R_3z, R_A, c1, c2, c3, c4, c5 = uks[0],uks[1],uks[2],uks[3],uks[4],uks[5],uks[6],uks[7],uks[8],uks[9],uks[10],uks[11]

#retriving functions for spanwise aerodunamic load for multiple integral terms in distribution expression
t_x    = ae.t_w
w_x    = ae.w_x
x_ae   = ae.x
z_ae   = ae.z
#print(x_ae)
#print(R_1y, R_2y, R_3y, R_1z, R_2z, R_3z, R_A, c1, c2, c3, c4, c5)

def mac(x, xi, n):
    r = (x-xi)      #this is the step function switch
    if r>0:
        return r**n
    else:
        return 0

def M_y(la, R_A, R_1z, R_2z, R_3z, x1, x2, x3, xaj, xp, theta, P):

    x = np.linspace(0,la,40)
    y = []

    for i in x:
        y.append(-R_1z*mac(i,x1,1)-R_2z*mac(i,x2,1)-R_3z*mac(i,x3,1) + P*m.cos(theta)*mac(i,xp,1) - R_A*m.cos(theta)*mac(i,xaj,1))

    return np.array(x), np.array(y)

def M_z(la, R_A, R_1y, R_2y, R_3y, x1, x2, x3, xaj,xp, theta, P):
    x = np.linspace(0, la, 40)
    y = []

    lift_object = ii.Interpolate_Integrate(x_ae, w_x)

    for i in x:
        y.append(-R_1y*mac(i,x1,1)  + lift_object.int_spline_natural(2,min(i,la_lim)) - R_2y*mac(i,x2,1)- R_3y*mac(i,x3,1) + P*m.sin(theta)*mac(i,xp,1)-R_A*m.sin(theta)*mac(i,xaj, 1))
        if i > x3:
            y[-1] = 0
    return np.array(x), np.array(y)


def S_y(la, R_A, R_1y, R_2y, R_3y, x1, x2, x3, xaj, xp, theta, P):
    x = np.linspace(0, la, 40)
    y = []

    lift_object = ii.Interpolate_Integrate(x_ae, w_x)

    for i in x:
        y.append(-R_1y*mac(i,x1,0)  + lift_object.int_spline_natural(1,min(i,la_lim)) - R_2y*mac(i,x2,0)- R_3y*mac(i,x3,0) + P*m.sin(theta)*mac(i,xp,0)-R_A*m.sin(theta)*mac(i,xaj, 0))

    return np.array(x), np.array(y)

def S_z(la, R_A, R_1z, R_2z, R_3z, x1, x2, x3, xaj, xp, theta, P):
    x = np.linspace(0, la, 40)
    y = []

    for i in x:
        y.append(-R_1z*mac(i,x1,0)-R_2z*mac(i,x2,0)-R_3z*mac(i,x3,0) + P*m.cos(theta)*mac(i,xp,0) - R_A*m.cos(theta)*mac(i,xaj,0))

    return np.array(x), np.array(y)

def v_x(la, R_A, R_1y, R_2y, R_3y, x1, x2, x3, xaj, xp, theta, P):
    x = np.linspace(0, la, 40)
    y = []

    lift_object = ii.Interpolate_Integrate(x_ae, w_x)

    for i in x:
        y.append((-1/E/I_zz)*(((-R_1y/6)*mac(i,x1,3)) + lift_object.int_spline_natural(4,min(i,la_lim/1.5)) - (R_2y/6)*mac(i,x2,3) -
                   (R_3y/6)*mac(i,x3,3) + (P/6)*np.sin(theta)*mac(i,xp,3) - (R_A/6)*np.sin(theta)*mac(i,xaj,3) ) + c1*i + c2)

    return np.array(x), np.array(y)

def other_x(la, R_A, R_1z, R_2z, R_3z, x1, x2, x3, xaj, xp, theta, P):  #this is  w(x)
    x = np.linspace(0, la, 40)
    y = []

    for i in x:
        y.append((-1/E/I_yy)*((-R_1z/6)*mac(i,x1,3) + (P/6)*np.cos(theta)*mac(i,xp,3)-(R_2z/6)*mac(i,x2,3)-(R_3z/6)*mac(i,x3,3) - (R_A/6)*np.cos(theta)*mac(i,xaj,3)) + c3*i + c4)

    return np.array(x), np.array(y)

def T_x(la, R_A, R_1y, R_2y, R_3y, x1, x2, x3, xaj, xp, theta, P):  #spanwise torque distribution

    x = np.linspace(0, la, 40)
    y = []

    torque_object = ii.Interpolate_Integrate(x_ae, t_x)

    for i in x:
        y.append(torque_object.int_spline_natural(1,min(i,la_lim)) + (ha/2 + z_sc)*R_1y*mac(i,x1,0) +
            (ha/2 + z_sc)*R_2y*mac(i,x2,0)+ (ha/2 + z_sc)*R_3y*mac(i,x3,0) + z_sc*R_A*m.sin(theta)*mac(i,xaj,0) +
            r * R_A * m.cos(theta) * mac(i, xaj, 0) -z_sc*P*m.sin(theta)*mac(i,xp,0)  - r*P*m.cos(theta) * mac(i, xp, 0))

    return np.array(x), np.array(y)


def th_x(la, R_A, R_1y, R_2y, R_3y, x1, x2, x3, xaj, xp, theta, P):  #spanwise torque distribution

    x = np.linspace(0, la, 40)
    y = []

    torque_object = ii.Interpolate_Integrate(x_ae, t_x)

    for i in x:
        #print(torque_object.int_spline_natural(2,min(i,la_lim/1.1)))
        if i < x3:
            y.append((1/G/J)*(torque_object.int_spline_natural(2,min(i,la_lim/1.1)) + (ha/2 + z_sc)*R_1y*mac(i,x1,1) +
                          (ha/2 + z_sc)*R_2y*mac(i,x2,1) + (ha/2 + z_sc)*R_3y*mac(i,x3,1) + z_sc*R_A*m.sin(theta)*mac(i,xaj,1) +
            r * R_A * m.cos(theta) * mac(i, xaj, 1) -z_sc*P*m.sin(theta)*mac(i,xp,1)  - r*P*m.cos(theta) * mac(i, xp, 1))+c5)
        else:
            y.append(y[-1])


    return np.array(x), np.array(y)


x_sy, y_sy = S_y(la, R_A, R_1y, R_2y, R_3y, x1, x2, x3, xaj, xp, theta, P)
x_sz, y_sz = S_z(la, R_A, R_1z, R_2z, R_3z, x1, x2, x3, xaj, xp, theta, P)
x_my, y_my = M_y(la, R_A, R_1z, R_2z, R_3z, x1, x2, x3, xaj, xp, theta, P)
x_mz, y_mz = M_z(la, R_A, R_1y, R_2y, R_3y, x1, x2, x3, xaj,xp, theta, P)
x_v,  y_v  = v_x(la, R_A, R_1y, R_2y, R_3y, x1, x2, x3, xaj, xp, theta, P)
x_w,  y_w  = other_x(la, R_A, R_1z, R_2z, R_3z, x1, x2, x3, xaj, xp, theta, P)
x_T,  y_T  = T_x(la, R_A, R_1y, R_2y, R_3y, x1, x2, x3, xaj, xp, theta, P)
x_th, y_th = th_x(la, R_A, R_1y, R_2y, R_3y, x1, x2, x3, xaj, xp, theta, P)


# print(x_ae)
# print(w_x)
# fig,a =  plt.subplots(2,3)
# a[0][0].plot(x_sy,y_sy)
# a[0][0].set_title('Shear in y')
# a[0][1].plot(x_sz,y_sz)
# a[0][1].set_title('Shear in z')
# a[0][2].plot(x_my,y_my)
# a[0][2].set_title('Moment in y')
# a[1][0].plot(x_mz,y_mz)
# a[1][0].set_title('Moment in z')
# a[1][1].plot(x_v,y_v)
# a[1][1].set_title('v(x)')
# a[1][2].plot(x_w,y_w)
# a[1][2].set_title('w(x)')
# plt.show()

# fig,b = plt.subplots(1,2)
# b[0][0].plot(x_T,y_T)
# b[0][0].set_title('Spanwise Torque Distribution')
# b[0][1].plot(x_th,y_th)
# b[0][1].set_title('Spanwise Twist Distribution')
# plt.plot(x_th,y_th)
# plt.show()
# plt.plot(x_T,y_T)
# plt.show()
