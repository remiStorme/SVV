import matplotlib.pyplot as plt
import numpy as np


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

        ii = np.digitize(self.xx, self.x_data)
        ii = np.fmin(np.fmax(ii - 1, 0), self.subspaces - 1)
        ff = np.zeros(self.xx.shape)

        for j, i in enumerate(ii):
            ff[j] = S(self.xx[j], self.x_data[i], a[i], b[i], c[i], d[i])

        return ff

    def plot_spline(self):
        yi = self.spline_natural()
        return plt.plot(self.xx, yi)
