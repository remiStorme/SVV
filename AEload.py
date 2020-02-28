import matplotlib.pyplot as plt
import numpy as np
import Interpolator_Integrate_Cubic as ii


plt.style.use('ggplot')

z_sc    = -0.08499063497059493
J       = 7.649955726444055*10**-6
I_yy    = 4.5925790464352304*10**-5
I_zz    = 4.686331664359035*10**-6



Nz = 81
Nx = 41
ca = 0.505
la = 1.611

#_____Fetching aerodynamic file

#including convertion to N/m**2
ae_data = 1000*np.genfromtxt("aerodynamicloadf100.dat", delimiter = ',')

#_____Creating angles for mesh coordinates
th_x  = np.zeros(Nx + 1)
th_z  = np.zeros(Nz + 1)

for i in range(Nx+1):
    th_x[i] = i*np.pi/Nx
for i in range(Nz+1):
    th_z[i] = i*np.pi/Nz


#_____Finding corrected mesh coordinates
x   = np.zeros(Nx)
z   = np.zeros(Nz)

for i in range(Nx):
    x[i]    = (la*(1-np.cos(th_x[i]))/2+la*(1-np.cos(th_x[i+1]))/2)/2

for i in range(Nz):
    z[i]    = (ca*(1-np.cos(th_z[i]))/2+ca*(1-np.cos(th_z[i+1]))/2)/2

#_____SPANWISE LIFT AND TORQUE DISTRIBUTION
#_____this is the functions that will be integrated in matrix b of the reaction forces system
#_____in order to construct these functions, the lift and torque due to it is integrated in z at
#each spanwise station
w_x = np.zeros(Nx)          #resultant lift distribution value for each spanwise location
t_w = np.zeros(Nx)          #the torque due to the ae-load at each spanwise location


for i in range(Nx):
    chord_lifts   = ae_data[:,i]   #select relevant data for the chosen spanwise location
    wtotal_x      = 0              #resultant contribution of lift for that spanwise location
    torque_x      = 0              #integral of the torque

    object        = ii.Interpolate_Integrate(z,chord_lifts) #creating a continuous function for chordwise lift
    wtotal_x     += object.int_spline_natural(1, ca)
    w_x[i]        = wtotal_x


    torque_function = -chord_lifts*(z + z_sc)      #creating an array representing the torque contribution
                                                  #of each chordwise lift value
    torque_object   = ii.Interpolate_Integrate(z,torque_function)
    torque_x       += torque_object.int_spline_natural(1,ca)
    t_w[i]          = torque_x

# plt.plot(x,t_w)
# plt.show()
#





