import numpy as np
from MOI import Iy, Iz, locbooms
from SC import gety, trapezoid, qbs
import Interpolator_Integrate_Cubic as ii

# Input Parameters
Ca = 0.505
h = 0.16
tsk = 0.0011
tsp = 0.0024
tst = 0.0012
hst = 0.013
wst = 0.017
nst = 11
Izz = Iz
Iyy = Iy
r = h / 2
l = np.sqrt(r ** 2 + (Ca - r) ** 2)
circ = np.pi * r / 2
alpha = np.arcsin(r / l)
T = 1
G = 1
Bi = wst * tst + hst * tst
s,_ = gety()
A_i = np.pi * r ** 2 / 2
A_ii = h * (Ca - r) / 2

q_int = []
s = [list(val) for val in s]
num_i=[]
for i in range(6):
    num_i.append(trapezoid(qbs[i],s[i])[-1]) # the last element of the list that "trapezoid" generates is the value of the integral
# print(num_i)

# TC is the matrix in terms of q0,1 and q0,2 for the compatibility equation and the total torque equation
TC = np.array([[((2 * circ)/(A_i * tsk)) + ((2 * r)/(A_i * tsp)) + ((2 * r)/(A_ii * tsk)), (-(2 * r)/(A_i * tsp)) - ((2 * r)/(A_ii * tsp)) - ((2 * l)/(A_ii * tsk))], [2 * A_i, 2 * A_ii]])

P = np.array([[((num_i[2] + num_i[3])/(A_ii * tsk)) + ((num_i[1] + num_i[4])/(A_ii * tsp)) - ((num_i[0] + num_i[5])/(A_i * tsk)) - ((num_i[1] + num_i[4])/(A_i * tsp))], [1]])
# print(P)

S = np.linalg.solve(TC, P)

q01_n = S[0]
q02_n = S[1]


alpha_twist = ((1)/(2 * A_i * G)) * (((num_i[0] + num_i[5])/(tsk)) + ((num_i[1] + num_i[4])/(tsp)) + q01_n[0] * (((s[0][-1] + s[5][-1])/(A_i * tsk)) + ((s[1][-1] + s[4][-1])/(A_i * tsp))) - q02_n[0] * ((s[1][-1] + s[4][-1])/(A_i * tsp)))

J = ((T)/(G * alpha_twist))

print(J)




"""
q_int = []
s = [list(val) for val in s]

for i in range(len(s)):
    q_object = ii.Interpolate_Integrate(s[i],qbs[i])
    # print(q_int)
    lim = s[i][-1]
    temp_int = q_object.int_spline_natural(1,lim)
    q_int.append(temp_int)

print(q_int)
"""
"""

# Boom Locations
cir = np.pi * r / 2
P = cir + l
s = 0
a = P / 6  # spacing
y_1 = []
y_3 = []
while s < P:
    if s < cir:
        y = r * np.sin(s / r)
        s += a
        y_1.append(y)
    elif cir < s < P:
        y = (P - s) * r / l
        s += a
        y_3.append(y)
print(y_1)
print(y_3)

B1 = Bi * sum(y_1)
B3 = Bi * sum(y_3)
#print(B1)
#print(B3)


y_1 = []
y_3 = []
for i in range(1):
    y_1.append(locbooms[i][1])
for j in range(6):
    y_3.append(locbooms[j][1])
B1 = Bi * sum(y_1)
B3 = Bi * sum(y_3)

# Tornsional stiffness
A_i = np.pi * r ** 2 / 2
A_ii = h * (Ca - r) / 2



D = np.array([[-1 * ((np.pi * r / t) + 2 * r / tsp) - 2 * r / t, -1 * (-2 * r / tsp) + l / t + 2 * r / tsp],
              [2 * A_i, 2 * A_ii]])

E = np.array([[-r ** 3 / (3 * Izz) + (1 / (Izz * t)) * (-2 * t * r ** 3 + (2 * B1) * r * np.pi / 2) + -1 * (
            r ** 3 / (3 * Izz) + (1 / (Izz * t)) * (
                2 * t * r * l ** 2 / 3 + 2 * t * r ** 2 * l + tsp * r ** 2 * l + 2 * (B1 + B3) * l + B3 * l))],
              [1]])

F = np.linalg.solve(D, E)

q01_n = F[0]
q02_n = F[1]

J = 1 / (G * (1 / (t * Izz) * (2 * t * r ** 3 - 2 * B1 * r * np.pi / 2) + r ** 3 / (3 * Izz) + q01_n * (
            np.pi * r / t + 2 * r / tsp) - 2 * r * q02_n))

print(J)
print(qbs)
"""