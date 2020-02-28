import numpy as np
from MOI import Iy, Iz
from SC import getszyf, trapezoid, qbs
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
s,_,_,_,_ = getszyf()
A_i = np.pi * r ** 2 / 2
A_ii = h * (Ca - r) / 2

q_int = []
s = [list(val) for val in s]

num_i= []
for i in range(6):
    num_i.append(trapezoid(qbs[i],s[i])[-1]) # the last element of the list that "trapezoid" generates is the value of
    # the integral

#Now we find the torsional stiffness

D = np.array([[-np.pi*r/(tsk * A_i) - 2*r/(tsp * A_i)- 2*r/(tsp * A_ii), 2*l/(tsk * A_ii) + 2*r/(tsp * A_ii)
               + 2*r/(tsp * A_i)],
              [2 * A_i, 2 * A_ii]])

E = np.array([[num_i[0]/tsk - 2*num_i[1]/tsp + 2*num_i[4]/tsp - num_i[5]/tsk - num_i[2]/tsk - num_i[3]/tsk],
           [T]])

F = np.linalg.solve(D, E)

q01n = F[0]
q02n = F[1]

twist1 = (1/2/A_i)*(num_i[0]/tsk - num_i[1]/tsp + num_i[4]/tsp - num_i[5]/tsk +(np.pi*r/tsk+2*r/tsp)*q01n - 2*r/tsp*q02n)
twist2 = (1/2/A_ii)*(num_i[1]/tsp + num_i[2]/tsk + num_i[3]/tsk - num_i[4]/tsp + (2*l/tsk+2*r/tsp)*q02n - 2*r/tsp*q01n)
torque = 2 * A_i * q01n + 2 * A_ii * q02n

#print(torque)
#print(twist1)
#print(twist2)


J = T / twist1
print('J =', J)

# Error compared to the answer of the verification model

e = ((J - 7.748548555816593 * 10**(-6))/(7.748548555816593 * 10**(-6))) * 100 # in percent
print('The error of J is =', e[0],'%')