import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, griddata, RegularGridInterpolator
from scipy.signal import savgol_filter
from scipy.constants import mu_0

#function definitions. 
def swap_order(X, I):
    #X is the array we want to swap the order of
    #I is the index of the element that we want to be first. 
    
    N = len(X)
    A = np.asarray(X[I:], float)
    B = np.asarray(X[0:I], float)
    Y = np.concatenate((A, B))
    return(Y)
#Your code may nneed to adjust the pressure gradient to get a solution, gives you the stability boundaries. 


#coil parameters

coil_height = 2.0
coil_major_radius = 1.0
coil_minor_radius = 0.1

#Upsampling and smoothing settings. 
N_dense = 10000
window_length = int((N_dense/40) + 1)*10  # Window size (must be odd)
polyorder = 2  # Polynomial order (usually 2 or 3)

# Loading the data -_-_-_-
pressure = np.loadtxt("pressure.txt")
psi = np.loadtxt("psi.txt")


pressure_max = 689476.0
psi_max = 1.0


R = np.loadtxt("R_grid.txt")
Z = np.loadtxt("Z_grid.txt")
NR = len(R)
NZ = len(Z)


R_grid, Z_grid = np.meshgrid(R, Z, indexing='ij')
DR = R[1]-R[0]
DZ = Z[1]-Z[0]

pi = np.pi

#Extracting Contours _-_-_-_

psi_level = (psi.max() - psi.min())*0.8#psi value of the contour I care about. 
d_psi = (psi.max() - psi.min())/300 #offset for grad psi. 
fig, ax = plt.subplots()
contour_1 = ax.contour(R_grid, Z_grid, psi, levels=[psi_level])
contour_path = contour_1.collections[0].get_paths()
contour_points = contour_path[0].vertices  # (N, 2) array of (R, Z) points
R_contour, Z_contour = contour_points[:, 0], contour_points[:, 1] #R, Z values of each point on the contour. 
R_indices = (R_contour/DR) #Needed to get the field points along the contour. 
Z_indices = (Z_contour/DZ)
N = len(R_indices)



#Checking to make sure the contour is fully periodic. 
print("R contour periodicity check: " + str((R_contour[-1] - R_contour[0])/(R_contour[0])))
print("Z contour periodicity check: " + str((Z_contour[-1] - Z_contour[0])/Z_contour[0]))



for n in range(0, N):
    R_indices[n] = int(R_indices[n])
    Z_indices[n] = int(Z_indices[n])


#calculating psi on the original less fine grid to later be interpolated on the contour.
d_psi_d_Z, d_psi_d_R = np.gradient(psi, DZ, DR)

d_psi_d_Z = savgol_filter(d_psi_d_Z, window_length, polyorder, mode='wrap')
d_psi_d_R = savgol_filter(d_psi_d_R, window_length, polyorder, mode='wrap')


d_psi_d_Z_contour = np.zeros(N)
d_psi_d_R_contour = np.zeros(N)

for n in range(0, N):
    R_ind_temp = int(R_indices[n])
    Z_ind_temp = int(Z_indices[n])
    d_psi_d_R_contour[n] = d_psi_d_R[R_ind_temp][Z_ind_temp]
    d_psi_d_Z_contour[n] = d_psi_d_Z[R_ind_temp][Z_ind_temp]

print("dpsi_dR contour periodicity check: " + str((d_psi_d_R_contour[-1] - d_psi_d_R_contour[0])/(d_psi_d_R_contour[0])))
print("dpsi_dZ contour periodicity check: " + str((d_psi_d_Z_contour[-1] - d_psi_d_Z_contour[0])/d_psi_d_Z_contour[0]))


#Calculating relevant quantities along the contour. -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

#l along the contour. 
l = np.zeros(N)
L = 0
for n in range(1, N):
    L += np.sqrt(((R_contour[n] - R_contour[n-1])**2.0) + ((Z_contour[n] - Z_contour[n-1])**2.0))
    l[n] = L

print("R contour periodicity check: " + str((R_contour[-1] - R_contour[0])/(R_contour[0])))
print("Z contour periodicity check: " + str((Z_contour[-1] - Z_contour[0])/Z_contour[0]))


spline_R = CubicSpline(l, R_contour, bc_type='periodic')
spline_Z = CubicSpline(l, Z_contour, bc_type='periodic')




#Defining the new denser points on the contour.
l_dense = np.linspace(0, L, N_dense)
R_dense = spline_R(l_dense)
Z_dense = spline_Z(l_dense)



#\grad \psi Calculation _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
#spline then filter, or filter then spline?
spline_dpsi_dR = CubicSpline(l, d_psi_d_R_contour, bc_type='periodic')
spline_dpsi_dZ = CubicSpline(l, d_psi_d_Z_contour, bc_type='periodic')
dpsi_dR_dense = savgol_filter(spline_dpsi_dR(l_dense), window_length, polyorder, mode='wrap')
dpsi_dZ_dense = savgol_filter(spline_dpsi_dZ(l_dense), window_length, polyorder, mode='wrap')
grad_psi = np.sqrt((dpsi_dR_dense**2.0) + (dpsi_dZ_dense**2.0))
grad_psi_2 = (dpsi_dR_dense**2.0) + (dpsi_dZ_dense**2.0)


print("dpsi_dR_dense contour periodicity check: " + str((dpsi_dR_dense[-1] - dpsi_dR_dense[0])/(dpsi_dR_dense[0])))
print("dpsi_dZ_dense contour periodicity check: " + str((dpsi_dZ_dense[-1] - dpsi_dZ_dense[0])/dpsi_dZ_dense[0]))

dpsi_dR_dense[-1] = dpsi_dR_dense[0]
dpsi_dZ_dense[-1] = dpsi_dZ_dense[0] #enforcing grad psi periodicity







#Magnetic Field Calculation _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
B_R = (-1.0)*dpsi_dZ_dense/R_dense
B_Z = dpsi_dR_dense/R_dense

#smoothing. 
B_R = savgol_filter(B_R, window_length, polyorder, mode='wrap')
B_R[0] = B_R[-1]
B_Z = savgol_filter(B_Z, window_length, polyorder, mode='wrap')
B_Z[0] = B_Z[-1]


B = np.sqrt((B_R**2.0) + (B_Z**2.0))
B_2 = (B_R**2.0) + (B_Z**2.0)




#checking to make sure that upsampling hasn't messed with periodicity. 
print("L dense periodicity check: " + str(L-l_dense[-1]))
print("R dense periodicity check: " + str((R_dense[-1] - R_dense[0])/(R_dense[0])))
print("Z dense periodicity check: " + str((Z_dense[-1] - Z_dense[0])/Z_dense[0]))


#Curvature Calculation _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- 

#Getting the derivatives that I need for curvature. 
R1_dense = spline_R.derivative(1)(l_dense)
Z1_dense = spline_Z.derivative(1)(l_dense)
R2_dense = spline_R.derivative(2)(l_dense)
Z2_dense = spline_Z.derivative(2)(l_dense)

print("pre smoothing")
print("R1 dense periodicity check: " + str((R1_dense[-1] - R1_dense[0])/(R1_dense[0])))
print("Z1 dense periodicity check: " + str((Z1_dense[-1] - Z1_dense[0])/Z1_dense[0]))

#smoothing the derivatives to avoid the sampling bumps. 
R1_dense = savgol_filter(R1_dense, window_length, polyorder, mode='wrap')
Z1_dense = savgol_filter(Z1_dense, window_length, polyorder, mode='wrap')
R2_dense = savgol_filter(R2_dense, window_length, polyorder, mode='wrap')
Z2_dense = savgol_filter(Z2_dense, window_length, polyorder, mode='wrap')


print("post smoothing")
print("R1 dense periodicity check: " + str((R1_dense[-1] - R1_dense[0])/(R1_dense[0])))
print("Z1 dense periodicity check: " + str((Z1_dense[-1] - Z1_dense[0])/Z1_dense[0]))


#enforcing periodicity
R1_dense[-1] = R1_dense[0]
R2_dense[-1] = R2_dense[0]

Z1_dense[-1] = Z1_dense[0]
Z2_dense[-1] = Z2_dense[0]

#tangent vector is \vec{R1_dense, Z1_dense}
unit_tangent = np.sqrt((R1_dense**2.0)+(Z1_dense**2.0))
normal_mag = np.sqrt((R2_dense**2.0)+(Z2_dense**2.0))
unit_normal_R = R2_dense/normal_mag
unit_normal_Z = Z2_dense/normal_mag

unit_normal_R = savgol_filter(unit_normal_R, window_length, polyorder, mode='wrap')
unit_normal_Z = savgol_filter(unit_normal_Z, window_length, polyorder, mode='wrap')
unit_normal_R[-1] = unit_normal_R[0]
unit_normal_Z[-1] = unit_normal_Z[0]



#curvature. _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
kappa_dense = ((R1_dense*Z2_dense) - (Z1_dense*R2_dense)) / ((R1_dense**2) + (Z1_dense**2))**1.5
kappa_dense_R = kappa_dense*unit_normal_R
kappa_dense_Z = kappa_dense*unit_normal_Z

kappa_dense = savgol_filter(kappa_dense, window_length, polyorder, mode='wrap')
kappa_dense_R = savgol_filter(kappa_dense_R, window_length, polyorder, mode='wrap')
kappa_dense_Z = savgol_filter(kappa_dense_Z, window_length, polyorder, mode='wrap')

kappa_dense[-1] = kappa_dense[0]
kappa_dense_R[-1] = kappa_dense_R[0]
kappa_dense_Z[-1] = kappa_dense_Z[0]



#grad psi dot kappa _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
grad_psi_dot_kappa = (kappa_dense_R*dpsi_dR_dense) + (kappa_dense_Z*dpsi_dR_dense)


closest_index = 0
closest_distance = 100
B_min = 100000
B_test = 100000
#Finding the spot closest to the surface of the ring:
for i_l in range(0, len(l_dense)):
    B_test = B[i_l]
    #d = np.sqrt(((coil_major_radius - R_dense[i_l])**2.0) + ((coil_height - Z_dense[i_l])**2.0))
    if B_test < B_min:
        B_min = B_test
        closest_index = i_l
    #if d < closest_distance:
        # closest_distance = d
        # closest_index = i_l
print("Min field " + str(B_min))
print("at index "+str(closest_index))
#Rearranging all data so that this closest index is the start, so that I can set ICs there and integrate along. 

R_dense = swap_order(R_dense, closest_index)
R1_dense = swap_order(R2_dense, closest_index)
R1_dense = swap_order(R2_dense, closest_index)
Z_dense = swap_order(Z_dense, closest_index)
Z1_dense = swap_order(Z1_dense, closest_index)
Z2_dense = swap_order(Z2_dense, closest_index)
B = swap_order(B, closest_index)
grad_psi = swap_order(grad_psi, closest_index)
grad_psi_2 = swap_order(grad_psi_2, closest_index)
kappa_dense = swap_order(kappa_dense, closest_index)
kappa_dense_R = swap_order(kappa_dense_R, closest_index)
kappa_dense_Z = swap_order(kappa_dense_Z, closest_index)
grad_psi_dot_kappa = swap_order(grad_psi_dot_kappa, closest_index)

def pressure_value(flux):
    return((pressure_max*(1 - np.cos((flux/psi_max)*np.pi))))
    
def pressure_derivative(flux):
    return((pressure_max*np.pi/psi_max)*(np.sin((flux)*np.pi/psi_max)))

dp_dpsi = pressure_derivative(psi_level)
print("pressure gradient: " + str(dp_dpsi))

#Now all the data is ready to be integrated. I start off assuming that the pertubation and it's first derivative
# is zero at the point closest to the coil, but that the second derivative has a very small positive value. 
starting_offset = 1e-9
Gamma = np.zeros(N_dense)

dl = L/N_dense
def A_1(l):
    return(grad_psi_2[l])

def A_2(l):
    return(2*grad_psi[l])

def A_3(l):
    return(((2*mu_0)/B_2[l])*(grad_psi_2[l])*(dp_dpsi)*(grad_psi_dot_kappa[l]))

def single_step(l):
    A1 = A_1(l)
    A2 = A_2(l)
    A3 = A_3(l)
    t1 = (((A2*dl)/2)-A1)
    t2 = ((2*A1)-(A3*dl*dl))
    t3 = (A1+((A2*dl)/2))

    return((t1+t2)/t3)
#starting the stepping off, we need forward differences. 
Gamma[1] = starting_offset
Gamma[3] = starting_offset

#forward stepping through the contour. 
for l in range(2, N_dense-1):
    Gamma[l+1] = single_step(l)


A_1_array = np.zeros(N_dense)
A_2_array = np.zeros(N_dense)
A_3_array = np.zeros(N_dense)

t_1_array = np.zeros(N_dense)
t_2_array = np.zeros(N_dense)
t_3_array = np.zeros(N_dense)

for l in range(0, N_dense):
    A_1_array[l] = A_1(l)
    A_2_array[l] = A_2(l)
    A_3_array[l] = A_3(l)
    A1 = A_1(l)
    A2 = A_2(l)
    A3 = A_3(l)
    t_1_array = (((A2*dl)/2)-A1)
    t_2_array = ((2*A1)-(A3*dl*dl))
    t_3_array = (A1+((A2*dl)/2))
    
    
A_1_array = savgol_filter(A_1_array, window_length, polyorder, mode='wrap')    
A_2_array = savgol_filter(A_2_array, window_length, polyorder, mode='wrap')    
A_3_array = savgol_filter(A_3_array, window_length, polyorder, mode='wrap')    



print("gamma = " + str(Gamma[100]))        
plt.figure()
plt.plot(l_dense, Gamma)
plt.title("Gamma")
plt.figure()
plt.plot(l_dense, A_1_array)
plt.title("A_1")
plt.figure()
# plt.plot(l_dense, A_2_array)
# plt.title("A_2")
# plt.figure()
# plt.plot(l_dense, A_3_array)
# plt.title("A_3")
# plt.figure()
plt.plot(l_dense, ((((A_2_array*dl)/2)-A_1_array)))
plt.title("t_1")
plt.figure()
plt.plot(l_dense, ((2*A_1_array)-(A_3_array*dl*dl)))
plt.title("t_2")
plt.figure()
plt.plot(l_dense, (A_1_array+((A_2_array*dl)/2)))
plt.title("t_3")
plt.figure()
#Plotting _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
# plt.plot(R_dense[0], Z_dense[0], 'xb')
# plt.plot(R_dense[int(N_dense/10)], Z_dense[int(N_dense/10)], 'xr')
# plt.figure()
# plt.plot(l_dense, R_dense)
# plt.title('R(l)')
# plt.figure()
# plt.plot(l_dense, Z_dense)
# plt.title('Z(l)')
# plt.figure()
# plt.plot(l_dense, R1_dense)
# plt.title('DR/DL')
# plt.figure()
# plt.plot(l_dense, R2_dense)
# plt.title('D2R/DL2')
# plt.figure()
# plt.plot(l_dense, Z1_dense)
# plt.title('DZ/DL')
# plt.figure()
# plt.plot(l_dense, Z2_dense)
# plt.title('D2Z/DL2')
# plt.figure()
# plt.plot(l_dense, kappa_dense)
# plt.title('curvature')
# plt.figure()
# plt.plot(l_dense, unit_tangent)
# plt.title("tangent check")
# plt.figure()
# plt.plot(l_dense, normal_mag)
# plt.title("normal check")
# plt.figure()
# plt.plot(l_dense, kappa_dense_R)
# plt.title("R curvature")
# plt.figure()
# plt.plot(l_dense, kappa_dense_Z)
# plt.title("Z curvature")
# plt.figure()
# plt.plot(l_dense, B)
# plt.title(r"B")
# plt.figure()
# plt.plot(l_dense, B_2)
# plt.title(r"B_2")
# plt.figure()
# plt.plot(l_dense, grad_psi)
# plt.title(r"\nabla \psi")
# plt.figure()



"""
#Making a denser psi grid _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
N_dense_psi = 1000
R_dense_psi = np.linspace(R_grid.min(), R_grid.max(), N_dense_psi)
Z_dense_psi = np.linspace(Z_grid.min(), Z_grid.max(), N_dense_psi)
DR_dense_psi = R_dense_psi[1] - R_dense_psi[0]
DZ_dense_psi = Z_dense_psi[1] - Z_dense_psi[0]
R_dense_psi, Z_dense_psi = np.meshgrid(R_dense_psi, Z_dense_psi)  # Create 2D dense grids

# Flatten the R_grid and Z_grid for use in griddata
points = np.array([R_grid.flatten(), Z_grid.flatten()]).T
values = psi.flatten()
psi_dense = griddata(points, values, (R_dense_psi, Z_dense_psi), method='cubic')  # Use 'linear', 'nearest', or 'cubic'
psi_dense = psi_dense.reshape(N_dense_psi, N_dense_psi)
plt.figure(figsize=(8, 6))
CS = plt.contour(R_dense_psi, Z_dense_psi, psi_dense, levels=50)#np.linspace(0, 0.5, 50))
plt.colorbar(label="Psi")
# plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Contour Plot of Psi")
"""








plt.show()

    



    
        


    
    





#Plotting Original Psi and Pressure_-_-_-

# #Plotting Psi
# plt.figure(figsize=(8, 6))
# CS = plt.contour(R_grid, Z_grid, psi, levels=50)#np.linspace(0, 0.5, 50))
# plt.colorbar(label="Psi")
# # plt.clabel(CS, inline=1, fontsize=10)
# plt.xlabel("R")
# plt.ylabel("Z")
# plt.title("Contour Plot of Psi")


# #Plotting Pressure
# plt.figure(figsize=(8, 6))
# CS = plt.contour(R_grid, Z_grid, pressure, levels=50)#np.linspace(0, 0.5, 50))
# plt.colorbar(label="Pressure")
# # plt.clabel(CS, inline=1, fontsize=10)
# plt.xlabel("R")
# plt.ylabel("Z")
# plt.title("Contour Plot of Pressure")
# plt.show()