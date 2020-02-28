import numpy as np
import math as m
import matplotlib.pyplot as plt
import AEload as ae
import Interpolator_Integrate_Cubic as ii

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

#HINGE AND ACTUATOR LOCATIONS: TEMPORARY VARIABLES
xaj     = x2 - xa/2
xp      = x2 + xa/2

#SECTION PROPERTIES
z_sc    = -0.08499063497059493
J       = 7.649955726444055*10**-6
I_yy    = 4.5925790464352304*10**-5
I_zz    = 4.686331664359035*10**-6


A_t = 0.002686478946739162
yc = -4.339264140786877e-19
zc = -0.24048766835061938


#MATERIAL PROPERTIES
E       = 72.9*10**9       # E-modulus (Pa)
G       = 27.1*10**9       # G-modulus (Pa)

#ADJUSTED MESH COORDINATED
x       = ae.x             #m
z       = ae.z             #m

#SPANWISE FUNCTIONS: LIFT ANF TORQUE
t_x    = ae.t_w
w_x    = ae.w_x

#creating an empty system Ax = b
#vector x = [Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Ar, C1, C2, C3, C4, C5]
A       = np.zeros((12,12))
b       = np.zeros((12,1))

#row 0 --> BC1 --> deflection(v) @ hinge1
A[0,0]  = 0
A[0,1]  = 0
A[0,2]  = 0
A[0,3]  = 0
A[0,4]  = 0
A[0,5]  = 0
A[0,6]  = 0
A[0,7]  = x1
A[0,8]  = 1
A[0,9]  = 0
A[0,10] = 0
A[0,11] = (ha/2 + z_sc)

#row 1 --> BC2 --> deflection(w) @ hinge1
A[1,0]  = 0
A[1,1]  = 0
A[1,2]  = 0
A[1,3]  = 0
A[1,4]  = 0
A[1,5]  = 0
A[1,6]  = 0
A[1,7]  = 0
A[1,8]  = 0
A[1,9]  = x1
A[1,10] = 1
A[1,11] = 0

#row 2 --> BC3 --> deflection(v) @ hinge2
A[2,0]  = (1/I_zz/E/6)*((x2-x1)**3)  + (1/G/J)*(ha/2 + z_sc)**2*(x2-x1)
A[2,1]  = 0
A[2,2]  = 0
A[2,3]  = 0
A[2,4]  = 0
A[2,5]  = 0
A[2,6]  = (1/E/I_zz/6)*m.sin(theta)*(x2-xaj)**3 + (ha/2 + z_sc)*(1/G/J)*(x2-xaj)*(z_sc*m.sin(theta) + r*m.cos(theta))
A[2,7]  = x2
A[2,8]  = 1
A[2,9]  = 0
A[2,10] = 0
A[2,11] = (ha/2 + z_sc)

#row 3 --> BC4 --> deflection(w) @ hinge2
A[3,0]  = 0
A[3,1]  = 0
A[3,2]  = 0
A[3,3]  = (1/E/I_yy/6)*(x2-x1)**3
A[3,4]  = 0
A[3,5]  = 0
A[3,6]  = (1/E/I_yy/6)*m.cos(theta)*(x2-xaj)**3
A[3,7]  = 0
A[3,8]  = 0
A[3,9]  = x2
A[3,10] = 1
A[3,11] = 0

#row 4 --> BC5 --> deflection(v) @ hinge3
A[4,0]  = (1/E/I_zz/6)*(x3-x1)**3 + ((ha/2 + z_sc)/G/J)*(x3-x1)
A[4,1]  = (1/E/I_zz/6)*(x3-x2)**3 + ((ha/2 + z_sc)/G/J)*(x3-x2)
A[4,2]  = 0
A[4,3]  = 0
A[4,4]  = 0
A[4,5]  = 0
A[4,6]  = (1/E/I_zz/6)*m.sin(theta)*(x3-xaj)**3 + 1/(G*J)*(ha/2 + z_sc)*(z_sc*m.sin(theta)*(x3-xaj)
                                                                         + r*m.cos(theta)*(x3-xaj))
A[4,7]  = x3
A[4,8]  = 1
A[4,9]  = 0
A[4,10] = 0
A[4,11] = (ha/2 + z_sc)

#row 5 --> BC6 --> deflection(w) @ hinge3
A[5,0]  = 0
A[5,1]  = 0
A[5,2]  = 0
A[5,3]  = (1/E/I_yy/6)*(x3-x1)**3
A[5,4]  = (1/E/I_yy/6)*(x3-x2)**3
A[5,5]  = 0
A[5,6]  = (1/E/I_yy/6)*m.cos(theta)*(x3-xaj)**3
A[5,7]  = 0
A[5,8]  = 0
A[5,9]  = x3
A[5,10] = 1
A[5,11] = 0

#row 6 --> BC7 --> deflection(w+c) @ jammed actuator xaj
A[6,0]  = (xaj-x1)**3*(m.sin(theta)/E/I_zz/6) + (xaj-x1)*(1/G/J)*(ha/2 + z_sc)*((ha/2)*m.cos(theta) + z_sc*m.sin(theta))
A[6,1]  = 0
A[6,2]  = 0
A[6,3]  = (xaj-x1)**3*(m.cos(theta)/E/I_yy/6)
A[6,4]  = 0
A[6,5]  = 0
A[6,6]  = 0
A[6,7]  = m.sin(theta)*xaj
A[6,8]  = m.sin(theta)
A[6,9]  = m.cos(theta)*xaj
A[6,10] = m.cos(theta)
A[6,11] = (ha/2)*m.cos(theta) + z_sc*m.sin(theta)

#row 7 --> BC_Fy
A[7,0]  = 1
A[7,1]  = 1
A[7,2]  = 1
A[7,3]  = 0
A[7,4]  = 0
A[7,5]  = 0
A[7,6]  = m.sin(theta)
A[7,7]  = 0
A[7,8]  = 0
A[7,9]  = 0
A[7,10] = 0
A[7,11] = 0

#row 8 --> BC_Fz
A[8,0]  = 0
A[8,1]  = 0
A[8,2]  = 0
A[8,3]  = 1
A[8,4]  = 1
A[8,5]  = 1
A[8,6]  = m.cos(theta)
A[8,7]  = 0
A[8,8]  = 0
A[8,9]  = 0
A[8,10] = 0
A[8,11] = 0

#row 9 --> BC_My
A[9,0]  = 0
A[9,1]  = 0
A[9,2]  = 0
A[9,3]  = -(la-x1)
A[9,4]  = -(la-x2)
A[9,5]  = -(la-x3)
A[9,6]  = -m.cos(theta)*(la-xaj)
A[9,7]  = 0
A[9,8]  = 0
A[9,9]  = 0
A[9,10] = 0
A[9,11] = 0

#row 10 --> BC_Mz
A[10,0]  = -(la-x1)
A[10,1]  = -(la-x2)
A[10,2]  = -(la-x3)
A[10,3]  = 0
A[10,4]  = 0
A[10,5]  = 0
A[10,6]  = -m.sin(theta)*(la-xaj)
A[10,7]  = 0
A[10,8]  = 0
A[10,9]  = 0
A[10,10] = 0
A[10,11] = 0

#row 11 --> BC_T(x)
A[11,0]  = (ha/2 + z_sc)
A[11,1]  = (ha/2 + z_sc)
A[11,2]  = (ha/2 + z_sc)
A[11,3]  = 0
A[11,4]  = 0
A[11,5]  = 0
A[11,6]  = z_sc*m.sin(theta) + r*m.cos(theta)
A[11,7]  = 0
A[11,8]  = 0
A[11,9]  = 0
A[11,10] = 0
A[11,11] = 0

#RHS
lift_object     = ii.Interpolate_Integrate(x,w_x)
torque_object   = ii.Interpolate_Integrate(x,t_x)

b[0,0]   = d1*np.cos(theta) + (1/I_zz/E)*(lift_object.int_spline_natural(4,x1)) - \
           (ha/2 + z_sc)*(1/G/J)*(torque_object.int_spline_natural(2,x1))

b[1,0]   = -d1*np.sin(theta)

b[2,0]   = (1/I_zz/E)*(lift_object.int_spline_natural(4,x2)) - (ha/2 + z_sc)*(1/G/J)\
           *(torque_object.int_spline_natural(2,x2))

b[3,0]   = 0

lift_object_temp = ii.Interpolate_Integrate(x, w_x * (x3 - x) ** 3 / 6)
torque_object_temp = ii.Interpolate_Integrate(x, t_x * (x3 - x) ** 1 / 2)
b[4,0]   = d3*np.cos(theta)+(1/I_zz/E)*((lift_object_temp.int_spline_natural(1,x3)) + (P/6)*np.sin(theta)*(x3 - xp)**3)\
           - (ha/2 + z_sc)*(1/G/J)*((torque_object_temp.int_spline_natural(2,min(x3,la_lim)))
                                    - z_sc*P*np.sin(theta)*(x3-xp) - r*P*np.cos(theta)*(x3-xp))

b[5,0]   = -d3*np.sin(theta) + (1/I_yy/E)*(P/6)*np.cos(theta)*(x3 - xp)**3

b[6,0]   = (np.sin(theta)/E/I_zz)*(lift_object.int_spline_natural(4,xaj)) \
           -torque_object.int_spline_natural(2,xaj)*(1/G/J)*((ha/2)*m.cos(theta) + z_sc*m.sin(theta))

b[7,0]   = lift_object.int_spline_natural(1,la) + P*np.sin(theta)

b[8,0]   = P*np.cos(theta)

b[9,0]   = -P*np.cos(theta)*(la-xp)

lift_object_temp = ii.Interpolate_Integrate(x, w_x * (la - x) ** 1 / 2)
b[10,0]  = -lift_object_temp.int_spline_natural(1,la) - P*np.sin(theta)*(la-xp)

b[11,0]  = -torque_object.int_spline_natural(1,la) + z_sc*P*np.sin(theta) + r*P*np.cos(theta)

uks = np.linalg.solve(A,b)
#print(reactions)
#print(lift_object.int_spline_natural(4,x1))
#print(lift_object.int_spline_natural(4,x2))
#print(lift_object.int_spline_natural(4,x3))

#vector x = [Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Ar, C1, C2, C3, C4, C5]

R_1y, R_2y, R_3y, R_1z, R_2z, R_3z, R_A, c1, c2, c3, c4, c5 = uks[0],uks[1],uks[2],uks[3],uks[4],uks[5],uks[6],uks[7],\
                                                              uks[8],uks[9],uks[10],uks[11]

print(R_1y, R_2y, R_3y, R_1z, R_2z, R_3z, R_A, c1, c2, c3, c4, c5)

#print(lift_object.int_spline_natural(2,la))

#print(uks)