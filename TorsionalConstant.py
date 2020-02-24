import numpy as np
from MOI import Iy, Iz, locbooms


# Input Parameters
Ca = 0.505
h = 0.16
t = 0.0011
tsp = 0.0024
tst = 0.0012
hst = 0.013
wst = 0.017
nst = 11
Izz = Iz
Iyy = Iy
r = h / 2
l = np.sqrt(r ** 2 + (Ca - r) ** 2)
alpha = np.arcsin(r / l)
T = 1  # to be changed once we have the interpolation
G = 1
Bi = wst * tst + hst * tst

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
"""

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