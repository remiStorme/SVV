import matplotlib.pyplot as plt
import numpy as np
from AEload import FileReader,Coordinates



class Interpolate:
    def __init__(self, x_data, y_data):
        self.x_data = x_data
        self.y_data = y_data
        self.x0 = self.x_data[0]
        self.xn = self.x_data[-1]
        self.delta_x = np.diff(x_data)
        self.delta_y = np.diff(y_data)
        self.size = len(x_data)
        self.subspaces = self.size - 1
        self.xx = np.linspace(self.x0, self.xn, 101)

    def sys_matrix(self):
        A = np.zeros((self.size, self.size))
        # boundary conditions
        A[0, 0] = 1
        A[-1, -1] = 1

        for i in range(1, self.subspaces):
            A[i, i] = 2 * (self.delta_x[i - 1] + self.delta_x[i])
            A[i, i - 1] = self.delta_x[i - 1]
            A[i, i + 1] = self.delta_x[i]

        return A

    def rhs_matrix(self):
        b = np.zeros(self.size)
        for i in range(1, self.subspaces):
            b[i] = 3. * (self.delta_y[i]) / self.delta_x[i] - 3 * (self.delta_y[i - 1]) / self.delta_x[i - 1]
        return b

    def M_coeffs(self):
        M = np.linalg.solve(self.sys_matrix(), self.rhs_matrix())
        return M

    def spline_coefficients(self):
        spline = np.zeros((self.subspaces, 4))
        M = self.M_coeffs()

        a = self.y_data
        b = (self.delta_y / self.delta_x) - self.delta_x / 3 * (2 * M[:-1] + M[1:])
        c = M
        d = (M[1:] - M[:-1]) / (3 * self.delta_x)

        return [a, b, c, d]

    def spline_natural(self):
        a, b, c, d = self.spline_coefficients()

        def S(x, xj, a, b, c, d):
            return a + b * (x - xj) + c * (x - xj) ** 2 + d * (x - xj) ** 3

        def IS(delta_x,a,b,c,d):
            return a * (delta_x) + (b/2) * (delta_x) ** 2 + (c/3) * (delta_x) ** 3 + (d/4) * (delta_x) ** 4

        ii = np.digitize(self.xx, self.x_data)
        ii = np.fmin(np.fmax(ii - 1, 0), self.subspaces - 1)
        ff = np.zeros(self.xx.shape)
        iff = np.zeros(self.xx.shape)
        ifff = np.zeros(self.xx.shape)

        for j, i in enumerate(ii):
            ff[j] = S(self.xx[j], self.x_data[i], a[i], b[i], c[i], d[i])
            ifff[j] = S(self.xx[j], self.x_data[i], a[i], b[i], c[i], d[i])*self.delta_x[i]
            iff[j] = IS(self.delta_x[i], a[i], b[i], c[i], d[i])



        return ff,sum(iff),sum(ifff)



    def plot_spline(self):
        yi,_ = self.spline_natural()
        fig = plt.figure()
        plt.plot(self.xx, yi,"r--")
        plt.show()


mat = FileReader("C:/Users/Paul Simon Schön/Downloads/aerodynamicloadf100.dat")
Z,X = Coordinates()


func = Interpolate(X,mat[:,0])
_,integral1,integral2 = func.spline_natural()

print(f"integral func*deltax:{integral1}\nintegral integrand:{integral2}")
