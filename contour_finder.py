import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, griddata, RegularGridInterpolator
from scipy.signal import savgol_filter



# Loading the data -_-_-_-
pressure = np.loadtxt("pressure.txt")
psi = np.loadtxt("psi.txt")

R = np.loadtxt("R_grid.txt")
Z = np.loadtxt("Z_grid.txt")
NR = len(R)
NZ = len(Z)


R_grid, Z_grid = np.meshgrid(R, Z, indexing='ij')
DR = R[1]-R[0]
DZ = Z[1]-Z[0]



#Extracting Contours _-_-_-_

psi_level = (psi.max() - psi.min())*0.8 #psi value of the contour I care about. 
d_psi = (psi.max() - psi.min())/300 #offset for grad psi. 
fig, ax = plt.subplots()
contour_1 = ax.contour(R_grid, Z_grid, psi, levels=[psi_level])
contour_path = contour_1.collections[0].get_paths()
contour_points = contour_path[0].vertices  # (N, 2) array of (R, Z) points
R_contour, Z_contour = contour_points[:, 0], contour_points[:, 1] #R, Z values of each point on the contour. 
R_indices = (R_contour/DR) #Needed to get the field points along the contour. 
Z_indices = (Z_contour/DZ)
N = len(R_indices)
plt.plot(R_contour[0], Z_contour[0], 'xb')
plt.figure()


#Checking to make sure the contour is fully periodic. 
print("R contour periodicity check: " + str((R_contour[-1] - R_contour[0])/(R_contour[0])))
print("Z contour periodicity check: " + str((Z_contour[-1] - Z_contour[0])/Z_contour[0]))



for n in range(0, N):
    R_indices[n] = int(R_indices[n])
    Z_indices[n] = int(Z_indices[n])


#calculating psi on the original less fine grid to later be interpolated on the contour.
d_psi_d_Z, d_psi_d_R = np.gradient(psi, DZ, DR)

d_psi_d_Z_contour = np.zeros(N)
d_psi_d_R_contour = np.zeros(N)

for n in range(0, N):
    R_ind_temp = int(R_indices[n])
    Z_ind_temp = int(Z_indices[n])
    d_psi_d_R_contour[n] = d_psi_d_R[R_ind_temp][Z_ind_temp]
    d_psi_d_Z_contour[n] = d_psi_d_Z[R_ind_temp][Z_ind_temp]



#Calculating relevant quantities along the contour. -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

#l along the contour. 
l = np.zeros(N)
L = 0
for n in range(1, N):
    L += np.sqrt(((R_contour[n] - R_contour[n-1])**2.0) + ((Z_contour[n] - Z_contour[n-1])**2.0))
    l[n] = L


#Upsampling and smoothing settings. 
N_dense = 100000
window_length = int((N_dense/30) + 1)  # Window size (must be odd)
polyorder = 3  # Polynomial order (usually 2 or 3)


spline_R = CubicSpline(l, R_contour, bc_type='periodic')
spline_Z = CubicSpline(l, Z_contour, bc_type='periodic')




#Defining the new denser points on the contour.
l_dense = np.linspace(0, L, N_dense)
R_dense = spline_R(l_dense)
Z_dense = spline_Z(l_dense)


#\grad \psi Calculation _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

spline_dpsi_dR = CubicSpline(l, d_psi_d_R_contour, bc_type='periodic')
spline_dpsi_dZ = CubicSpline(l, d_psi_d_Z_contour, bc_type='periodic')
dpsi_dR_dense = savgol_filter(spline_dpsi_dR(l_dense), window_length, polyorder, mode='wrap')
dpsi_dZ_dense = savgol_filter(spline_dpsi_dZ(l_dense), window_length, polyorder, mode='wrap')
grad_psi = np.sqrt((dpsi_dR_dense**2.0) + (dpsi_dZ_dense**2.0))


plt.plot(l_dense, grad_psi)
plt.title(r"\nabla \psi")
plt.figure()



#Magnetic Field Calculation _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
B_R = (-1.0)*dpsi_dZ_dense/R_dense
B_Z = dpsi_dR_dense/R_dense
B = np.sqrt((B_R**2.0) + (B_Z**2.0))

plt.plot(l_dense, B)
plt.title(r"B")
plt.figure()


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

#smoothing the derivatives to avoid the sampling bumps. 
R1_dense = savgol_filter(R1_dense, window_length, polyorder, mode='wrap')
Z1_dense = savgol_filter(Z1_dense, window_length, polyorder, mode='wrap')
R2_dense = savgol_filter(R2_dense, window_length, polyorder, mode='wrap')
Z2_dense = savgol_filter(Z2_dense, window_length, polyorder, mode='wrap')

#tangent vector is \vec{R1_dense, Z1_dense}
unit_tangent = np.sqrt((R1_dense**2.0)+(Z1_dense**2.0))
normal_mag = np.sqrt((R2_dense**2.0)+(Z2_dense**2.0))
unit_normal_R = R2_dense/normal_mag
unit_normal_Z = Z2_dense/normal_mag


#curvature. _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
kappa_dense = ((R1_dense*Z2_dense) - (Z1_dense*R2_dense)) / ((R1_dense**2) + (Z1_dense**2))**1.5
kappa_dense_R = kappa_dense*unit_normal_R
kappa_dense_Z = kappa_dense*unit_normal_Z



#Plotting _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
plt.plot(l_dense, R_dense)
plt.title('R(l)')
plt.figure()
plt.plot(l_dense, Z_dense)
plt.title('Z(l)')
plt.figure()
plt.plot(l_dense, R1_dense)
plt.title('DR/DL')
plt.figure()
plt.plot(l_dense, R2_dense)
plt.title('D2R/DL2')
plt.figure()
plt.plot(l_dense, Z1_dense)
plt.title('DZ/DL')
plt.figure()
plt.plot(l_dense, Z2_dense)
plt.title('D2Z/DL2')
plt.figure()
plt.plot(l_dense, kappa_dense)
plt.title('curvature')
plt.figure()
plt.plot(l_dense, unit_tangent)
plt.title("tangent check")
plt.figure()
plt.plot(l_dense, normal_mag)
plt.title("normal check")
plt.figure()
plt.plot(l_dense, kappa_dense_R)
plt.title("r curvature")
plt.figure()
plt.plot(l_dense, kappa_dense_Z)
plt.title("X curvature")
plt.figure()



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