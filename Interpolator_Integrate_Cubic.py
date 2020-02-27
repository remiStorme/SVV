import numpy as np
# import matplotlib.pyplot as plt

class Interpolate_Integrate:
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

    def int_spline_natural(self, integrals, x):
        a, b, c, d = self.spline_coefficients()

        def S(x, xj, a, b, c, d):
            return a + b*(x - xj) + c*(x - xj)** 2 + d*(x - xj)** 3


        def int_1(x, xi, a, b, c, d):
            return a*(x - xi) + (b / 2)*(x - xi)**2 + (c / 3)*(x - xi)** 3 + (d / 4)*(x - xi)** 4

        def int_2(x, xi, x0, a, b, c, d):
            C1 = a * (x0-xi) + (b / 2) * (x0 - xi) ** 2 + (c / 3) * (x0 - xi) ** 3 + (d / 4) * (x0 - xi) ** 4

            return (a / 2)*(x - xi)** 2 + (b / 6)*(x - xi)**3 + (c / 12)*(x - xi)** 4 + (d / 20)*(x - xi)**5 - C1*(x - xi)

        def int_3(x, xi, x0, a, b, c, d):
            C1 = a * (x0 - xi) + (b / 2) * (x0 - xi) ** 2 + (c / 3) * (x0 - xi) ** 3 + (d / 4) * (x0 - xi) ** 4
            C2 = (a / 2) * (x0-xi) ** 2 + (b / 6) * (x0 - xi) ** 3 + (c / 12) * (x0 - xi) ** 4 + (d / 20) * (x0 - xi) ** 5 - C1 * (x0-xi)

            return (a / 6)*(x - xi)**3 + (b / 24)*(x - xi)**4 + (c / 60)*(x - xi)**5 + (d / 120)*(x - xi)**6 - (C1/2)*(x - xi)**2 - C2*(x-xi)

        def int_4(x, xi, x0, a, b, c, d):
            C1 = a * (x0 - xi) + (b / 2) * (x0 - xi) ** 2 + (c / 3) * (x0 - xi) ** 3 + (d / 4) * (x0 - xi) ** 4
            C2 = (a / 2) * (x0-xi) ** 2 + (b / 6) * (x0 - xi) ** 3 + (c / 12) * (x0 - xi) ** 4 + (d / 20) * (x0 - xi) ** 5 - C1 * (x0-xi)
            C3 = (a / 6)*(x0-xi)**3 + (b/24)*(x0-xi)**4 + (c/60)*(x0-xi)**5 + (d/120)*(x0-xi)**6 - (C1/2)*(x0-xi)**2 - C2*(x0-xi)

            return (a/24)*(x-xi)**4 + (b/120)*(x-xi)**5 + (c/360)*(x-xi)**6 + (d/840)*(x-xi)**7 - (C1/6)*(x-xi)**3 - (C2/2)*(x-xi)**2 - C3*(x-xi)

        ii = np.digitize(x, self.x_data)
        ii = np.fmin(np.fmax(ii - 1, 0), self.subspaces - 1)


        int_sums = 0


        for i in range(ii):
            if integrals == 1:
                int_sums += int_1(self.x_data[i + 1], self.x_data[i], a[i], b[i], c[i], d[i])



            elif integrals == 2:
                int_sums += int_2(self.x_data[i + 1], self.x_data[i],self.x_data[0], a[i], b[i], c[i], d[i])


            elif integrals == 3:
                int_sums += int_3(self.x_data[i + 1], self.x_data[i],self.x_data[0], a[i], b[i], c[i], d[i])

            elif integrals == 4:
                int_sums += int_4(self.x_data[i + 1], self.x_data[i],self.x_data[0], a[i], b[i], c[i], d[i])


        if integrals == 1:
            int_sums += int_1(x, self.x_data[ii],a[ii], b[ii], c[ii], d[ii])


        elif integrals == 2:
            int_sums += int_2(x, self.x_data[ii],self.x_data[0], a[ii], b[ii], c[ii], d[ii])


        elif integrals == 3:
            int_sums += int_3(x, self.x_data[ii],self.x_data[0], a[ii], b[ii], c[ii], d[ii])

        elif integrals == 4:
            int_sums += int_4(x, self.x_data[ii],self.x_data[0], a[ii], b[ii], c[ii], d[ii])

        return int_sums

    # def plot_spline(self):
    #     yi = self.spline_natural()
    #     plt.plot(self.xx, yi)
    #     plt.show()

#
# x0,xn = -1,1.0
# subspaces = 40
# x_i = np.linspace(x0,xn,subspaces+1)
#
#
# #test function
# def f_runge(x):
#     return x**2
#
# y_i = f_runge(x_i)
#
#
# s = Interpolate_Integrate(x_i,y_i)
# s.plot_spline()
#
# int= s.int_spline_natural(2, 0.99)
# plt.plot(x_i[2:],int_parts)
# plt.show()
#
# print(int_parts)
