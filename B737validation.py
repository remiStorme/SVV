import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

#Functions to use 

def location_elements(fname):
        
    list_elements = []
    #Get the average location of each element
    list_average_elements = []


    location = pd.read_csv(fname, delimiter = ',', skiprows = 9, nrows = 6588, names = ['Node','x','y','z'])
    elements_list = pd.read_csv(fname, delimiter = ',', skiprows = 6598, nrows = 6633,
                                names = ['Element','N1','N2','N3','N4'])
    
    x = location['x'].values
    y = location['y'].values
    z = location['z'].values
    nodes = location['Node'].values
    Element = elements_list['Element'].values
    
    N1 = elements_list['N1'].values #List of the 1st entry of nodes
    N2 = elements_list['N2'].values #List of the 2th entry of nodes
    N3 = elements_list['N3'].values #List of the 3th entry of nodes
    N4 = elements_list['N4'].values #List of the 4th entry of nodes
    
    #Make the elements
    for i in range(len(Element)):
    
        entry1 = np.transpose(np.array([x[np.where(nodes == N1[i])], y[np.where(nodes == N1[i])],
                                        z[np.where(nodes == N1[i])]])) #XYZ locations of the nodes, so weget a 4x3 matrix
        entry2 = np.transpose(np.array([x[np.where(nodes == N2[i])], y[np.where(nodes == N2[i])],
                                        z[np.where(nodes == N2[i])]]))
        entry3 = np.transpose(np.array([x[np.where(nodes == N3[i])], y[np.where(nodes == N3[i])],
                                        z[np.where(nodes == N4[i])]]))
        entry4 = np.transpose(np.array([x[np.where(nodes == N4[i])], y[np.where(nodes == N3[i])],
                                        z[np.where(nodes == N4[i])]]))
    
        element = np.vstack((entry1,entry2,entry3,entry4))
        list_elements.append(element) #List of all the elements
    
    #Get the average location for each element
    for i in range(1,len(list_elements)-1):
        x_avg = sum(list_elements[i][:,0])/4
        y_avg = sum(list_elements[i][:,1])/4
        z_avg = sum(list_elements[i][:,2])/4
        average_location= [x_avg,y_avg,z_avg,i]
        list_average_elements.append(average_location) #Gives the average for the location of the node
        # (kind of like the centroid)

    list_average_elements = np.array(list_average_elements)
    sorted_avg_elements = list_average_elements[np.argsort(list_average_elements[:, 0])] #sorted list of elements
    # in the x-axis
    
    return sorted_avg_elements

def stresses(fname, start, end):
    
    stress = pd.read_csv(fname, delimiter = ',', skiprows = start, nrows = end,
                         names = ['Element','Int','s_mis1','s_mis2','s_shear1','s_shear2'])
    
    s_mis1 = stress['s_mis1'].values
    s_mis2 = stress['s_mis2'].values
    s_shear1 = stress['s_shear1'].values
    s_shear2 = stress['s_shear2'].values
    
    s_mis_mean = np.array((s_mis1+s_mis2)/2)
    s_shear_mean = np.array((s_shear1+s_shear2)/2)
    
    return s_mis_mean, s_shear_mean, stress

def full_element(fname_rpt,fname_inp,start,end):
    
    s_mis_mean, s_shear_mean, stress = stresses(fname_rpt,start,end)
    sorted_avg_elements = location_elements(fname_inp)
    
    uncomplete_list_elements = stress['Element'].values
    element_T = np.array(uncomplete_list_elements)
    
    #Get a nice matrix so that putting them at the right spot works
    stresses_matrix = np.transpose(np.vstack((s_mis_mean,s_shear_mean,element_T)))
    zeros = np.zeros((6631,2))
    sorted_avg_elements = np.append(sorted_avg_elements,zeros,axis = 1)
    
    #Placing the stress at the right element, since the stresses are already sorted by element 
    for i in range(0,len(stresses_matrix)):
        for j in range(0,len(sorted_avg_elements)):
            if sorted_avg_elements[j][3] == stresses_matrix[i][2]:
               
                sorted_avg_elements[j][4] = stresses_matrix[i][0] 
                sorted_avg_elements[j][5] = stresses_matrix[i][1]
    
    x_avg = sorted_avg_elements[:,0]
    y_avg = sorted_avg_elements[:,1]
    z_avg = sorted_avg_elements[:,2]
    
    vonmissen_mean = sorted_avg_elements[:,4]
    shear_mean = sorted_avg_elements[:,5]
    
    return x_avg, y_avg, z_avg, vonmissen_mean, shear_mean


def maximum(coordinate,stress):
    points = sorted(list(set(coordinate))) #Sorted list of unique values in x, y or z
    entries = []
    stresses = []
    stress_max = []
    
    for i in  range(len(points)):
        entries.append(np.where(coordinate == points[i])) #append points when value is equal to the unique value,
        # gives list of entries where value is the same
    for i in range(len(entries)):
        stresses.append(stress[entries[i]]) #Get a list of the stresses when the entries are the same
    for i in range(len(stresses)):
        stress_max.append(max(stresses[i])) #Find the maximum along that chord
    
    return points, stress_max

def coordinates(fname_rpt,fname_inp,start,end,coordinate,stress): #coordinate and stress are from 0 to 2
    # and 3 to 4 respecitivly
     element_values = full_element(fname_rpt,fname_inp,start,end)
     coordinate = element_values[coordinate] #x = 0, y = 1, z = 2
     stress = element_values[stress] #vonmissen = 3, shear = 4
     
     points,stress_max = maximum(coordinate,stress)
     
     return points, stress_max
 
#File names

fname_rpt = 'B737_report.csv'
fname_inp = 'B737_input.csv'
start = 6687
end = 5777


deflection = pd.read_csv(fname_rpt, delimiter = ',', skiprows = 26706, nrows = 6587,
                         names = ['Node','Magnitude','U1','U2','U3'])

x,stress = coordinates(fname_rpt,fname_inp,start,end,0,3)
z = coordinates(fname_rpt,fname_inp,start,end,2,3)[0]

plt.title('Stresses along the span of the wing')
plt.xlabel('x (mm)')
plt.ylabel('stress (Pa)')
plt.scatter(x,stress, label = 'Vonmissen')
plt.legend()

plt.show()









