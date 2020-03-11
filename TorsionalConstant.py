import numpy as np
from MOI import Iy, Iz, locbooms
from SC import getszyf, trapezoid, Ygetqbs
# import Interpolator_Integrate_Cubic as ii

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
s,_,_,_,_= getszyf()
A_i = np.pi * r ** 2 / 2
A_ii = h * (Ca - r) / 2

# q_int = []



#Now we find the torsional stiffness

def Tgetq(T):
    qbs = Ygetqbs(1)
    num_i = []
    for i in range(6):
        num_i.append(trapezoid(qbs[i], s[i])[-1])  # the last element of the list that "trapezoid" generates is the value of the integral
    D = np.array([[-np.pi*r/(tsk * A_i) - 2*r/(tsp * A_i)- 2*r/(tsp * A_ii), 2*l/(tsk * A_ii) + 2*r/(tsp * A_ii) + 2*r/(tsp * A_i)],
                  [2 * A_i, 2 * A_ii]])

    E = np.array([[num_i[0]/tsk - 2*num_i[1]/tsp + 2*num_i[4]/tsp - num_i[5]/tsk - num_i[2]/tsk - num_i[3]/tsk],
               [T]])

    F = np.linalg.solve(D, E)

    q01n = F[1][0]
    q02n = F[1][0]

    Tq=[]
    for i in range(len(s)):
        Tq_i=[]
        if i==0 or i==5:
            for j in range(len(qbs[i])):
                Tq_i.append(qbs[i][j]+q01n)
        if i==1:
            for j in range(len(qbs[i])):
                Tq_i.append(qbs[i][j]-q01n+q02n)
        if i==2 or i==3:
            for j in range(len(qbs[i])):
                Tq_i.append(qbs[i][j]+q02n)
        if i==4:
            for j in range(len(qbs[i])):
                Tq_i.append(qbs[i][j]+q01n-q02n)
        Tq.append(Tq_i)

        return num_i,q01n,q02n,Tq


num_i,q01n,q02n,Tq=Tgetq(T)


twist1 = (1/2/A_i)*(num_i[0]/tsk - num_i[1]/tsp + num_i[4]/tsp - num_i[5]/tsk +(np.pi*r/tsk+2*r/tsp)*q01n - 2*r/tsp*q02n)
twist2 = (1/2/A_ii)*(num_i[1]/tsp + num_i[2]/tsk + num_i[3]/tsk - num_i[4]/tsp + (2*l/tsk+2*r/tsp)*q02n - 2*r/tsp*q01n)

J = T / twist1
print(J)

# Error compared to the answer of the verification model

e = ((J - 7.649955726444055 * 10**(-6))/(7.649955726444055 * 10**(-6))) * 100 # in percent
# print('The error is =', e[0],'%')
'''
q_int = []
s = [list(val) for val in s]
num_i= []
for i in range(6):
    num_i.append(trapezoid(qbs[i],s[i])[-1]) # the last element of the list that "trapezoid" generates is the value of the integral
# print(num_i)

D = np.array([[-1 * ((np.pi * r / tsk) + 2 * r / tsp) - 2 * r / tsk, -1 * (-2 * r / tsp) + l / tsk + 2 * r / tsp], [2 * A_i, 2 * A_ii]])

E = np.array([[- num_i[0] / tsk + num_i[1] / tsp - num_i[4] / tsp + num_i[5] / tsk + (num_i[2] / tsk + num_i[3] / tsk - num_i[4] / tsp + num_i[1] / tsp)], [T]])

F = np.linalg.solve(D, E)

q01n = F[0]
q02n = F[1]

J = T / ((1 / (2 * A_i))(num_i[0] / tsk - num_i[1] / tsp + num_i[4] / tsp + num_i[5] / tsk + q01n[0] * ((np.pi * r / tsk) + 2 * r / tsp) + q02n[0] * (-2 * r / tsp)))

print(J)
'''
'''
# TC is the matrix in terms of q0,1 and q0,2 for the compatibility equation and the total torque equation
TC = np.array([[((2 * circ)/(A_i * tsk)) + ((2 * r)/(A_i * tsp)) + ((2 * r)/(A_ii * tsk)), (-(2 * r)/(A_i * tsp)) - ((2 * r)/(A_ii * tsp)) - ((2 * l)/(A_ii * tsk))], [2 * A_i, 2 * A_ii]])

P = np.array([[((num_i[2] + num_i[3])/(A_ii * tsk)) + ((num_i[1] + num_i[4])/(A_ii * tsp)) - ((num_i[0] + num_i[5])/(A_i * tsk)) - ((num_i[1] + num_i[4])/(A_i * tsp))], [1]])
# print(P)

S = np.linalg.solve(TC, P)

q01_n = S[0]
q02_n = S[1]

alpha_twist = ((1)/(2 * A_i * G * 10**3)) * (((num_i[0] + num_i[5])/(tsk)) + ((num_i[1] + num_i[4])/(tsp)) + q01_n[0] * (((s[0][-1] + s[5][-1])/(A_i * tsk)) + ((s[1][-1] + s[4][-1])/(A_i * tsp))) - q02_n[0] * ((s[1][-1] + s[4][-1])/(A_i * tsp)))

J = ((T)/(G * alpha_twist))

print("J Method I =", J)

'''
'''
J1 = ((r * (np.pi/2) * tsk**3)/(3))
J2 = ((r * tsp**3)/(3))
J3 = ((l * tsk**3)/(3))
J4 = J3
J5 = J2
J6 = J1

JJ = J1 + J2 + J3 + J4 + J5 + J6
print("J Method II =", JJ)
'''


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
"""