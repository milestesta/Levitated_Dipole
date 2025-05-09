import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, griddata, RegularGridInterpolator
from scipy.signal import savgol_filter
from scipy.constants import mu_0
from scipy.integrate import solve_ivp

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

coil_height = 1.5
coil_major_radius = 0.4
coil_minor_radius = 0.1

#Upsampling and smoothing settings. 
N_dense = 10000
window_length = int((N_dense/10) + 1)*1  # Window size (must be odd)
polyorder = 2  # Polynomial order (usually 2 or 3)

# Loading the data -_-_-_-
pressure = np.loadtxt("pressure.txt")
psi = np.loadtxt("psi.txt")


pressure_max = 689476.0
psi_max = 0.08


R = np.loadtxt("R_grid.txt")
Z = np.loadtxt("Z_grid.txt")
NR = len(R)
NZ = len(Z)


R_grid, Z_grid = np.meshgrid(R, Z, indexing='ij')
DR = R[1]-R[0]
DZ = Z[1]-Z[0]

pi = np.pi

#Extracting Contours _-_-_-_

psi_level = (psi.max() - psi.min())*0.3#psi value of the contour I care about. 
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


# spline_R = CubicSpline(l, R_contour, bc_type='periodic')
# spline_Z = CubicSpline(l, Z_contour, bc_type='periodic')

space_window = N/30
space_poly = 2

# plt.figure()
# plt.plot(R_contour, Z_contour)
# plt.title('before')
R_contour = savgol_filter(R_contour, space_window, space_poly, mode='wrap')
Z_contour = savgol_filter(Z_contour, space_window, space_poly, mode='wrap')

R_contour[-1] = R_contour[0]
Z_contour[-1] = Z_contour[0]
spline_R = CubicSpline(l, R_contour, bc_type='periodic')
spline_Z = CubicSpline(l, Z_contour, bc_type='periodic')
# plt.figure()
# plt.plot(R_contour, Z_contour)
# plt.title('after')


#Defining the new denser points on the contour.

dense_space_window = N_dense/15
dense_space_poly = 2
l_dense = np.linspace(0, L, N_dense)
R_dense = spline_R(l_dense)
Z_dense = spline_Z(l_dense)
# plt.figure()
# plt.plot(R_dense, Z_dense)
# plt.title('before')

R_dense = savgol_filter(R_dense, dense_space_window, dense_space_poly, mode='wrap')
Z_dense = savgol_filter(Z_dense, dense_space_window, dense_space_poly, mode='wrap')
R_dense = savgol_filter(R_dense, dense_space_window, dense_space_poly, mode='wrap')
Z_dense = savgol_filter(Z_dense, dense_space_window, dense_space_poly, mode='wrap')
R_dense = savgol_filter(R_dense, dense_space_window, dense_space_poly, mode='wrap')
Z_dense = savgol_filter(Z_dense, dense_space_window, dense_space_poly, mode='wrap')
R_dense = savgol_filter(R_dense, dense_space_window, dense_space_poly, mode='wrap')
Z_dense = savgol_filter(Z_dense, dense_space_window, dense_space_poly, mode='wrap')
R_dense[-1] = R_dense[0]
Z_dense[-1] = Z_dense[0]

# plt.figure()
# plt.plot(R_dense, Z_dense)
# plt.title('after')


#\grad \psi Calculations _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
#spline then filter, or filter then spline?

dpsi_window = N/10
dpsi_order = 2
# plt.figure()
# plt.plot(d_psi_d_R_contour)
# plt.title("before")

d_psi_d_R_contour = savgol_filter(d_psi_d_R_contour, dpsi_window, dpsi_order, mode='wrap')
d_psi_d_Z_contour = savgol_filter(d_psi_d_Z_contour, dpsi_window, dpsi_order, mode='wrap')
d_psi_d_R_contour = savgol_filter(d_psi_d_R_contour, dpsi_window, dpsi_order, mode='wrap')
d_psi_d_Z_contour = savgol_filter(d_psi_d_Z_contour, dpsi_window, dpsi_order, mode='wrap')
d_psi_d_R_contour[-1] = d_psi_d_R_contour[0]
d_psi_d_Z_contour[-1] = d_psi_d_Z_contour[0]
# plt.figure()
# plt.plot(d_psi_d_R_contour)
# plt.title("after")
spline_dpsi_dR = CubicSpline(l, d_psi_d_R_contour, bc_type='periodic')
spline_dpsi_dZ = CubicSpline(l, d_psi_d_Z_contour, bc_type='periodic')
dpsi_dR_dense = spline_dpsi_dR(l_dense)
dpsi_dZ_dense = spline_dpsi_dZ(l_dense)
# plt.figure()
# plt.plot(dpsi_dR_dense)
# plt.title("before")

dpsi2_window = N_dense/4
dpsi2_order = 3

dpsi_dR_dense = savgol_filter(dpsi_dR_dense, dpsi2_window, dpsi2_order, mode='wrap')
dpsi_dZ_dense = savgol_filter(dpsi_dZ_dense, dpsi2_window, dpsi2_order, mode='wrap')
dpsi_dR_dense = savgol_filter(dpsi_dR_dense, dpsi2_window, dpsi2_order, mode='wrap')
dpsi_dZ_dense = savgol_filter(dpsi_dZ_dense, dpsi2_window, dpsi2_order, mode='wrap')

dpsi_dR_dense[-1] = dpsi_dR_dense[0]
dpsi_dZ_dense[-1] = dpsi_dZ_dense[0] #enforcing grad psi periodicity


# plt.figure()
# plt.plot(dpsi_dR_dense)
# plt.title("after")

grad_psi_window = N_dense/10
grad_psi_order = 2
grad_psi = np.sqrt((dpsi_dR_dense**2.0) + (dpsi_dZ_dense**2.0))
# plt.figure()
# plt.plot(grad_psi)
# plt.title("before")

grad_psi = savgol_filter(grad_psi, grad_psi_window, grad_psi_order, mode='wrap')
grad_psi = savgol_filter(grad_psi, grad_psi_window, grad_psi_order, mode='wrap')
grad_psi = savgol_filter(grad_psi, grad_psi_window, grad_psi_order, mode='wrap')
grad_psi[-1] = grad_psi[0]
# plt.figure()
# plt.plot(grad_psi)
# plt.title("after")
spline_grad_psi = CubicSpline(l_dense, grad_psi, bc_type='periodic') #so I can use .derivative. 
grad_psi = spline_grad_psi(l_dense)
grad_psi_2 = (dpsi_dR_dense**2.0) + (dpsi_dZ_dense**2.0)


print("dpsi_dR_dense contour periodicity check: " + str((dpsi_dR_dense[-1] - dpsi_dR_dense[0])/(dpsi_dR_dense[0])))
print("dpsi_dZ_dense contour periodicity check: " + str((dpsi_dZ_dense[-1] - dpsi_dZ_dense[0])/dpsi_dZ_dense[0]))




#Magnetic Field Calculation _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
B_R = (-1.0)*dpsi_dZ_dense/R_dense
B_Z = dpsi_dR_dense/R_dense

B_window = N_dense/10
B_order = 2
#smoothing. 
# plt.figure()
# plt.plot(B_R)
# plt.title("before")
B_R = savgol_filter(B_R, B_window, B_order, mode='wrap')
B_Z = savgol_filter(B_Z, B_window, B_order, mode='wrap')
B_R = savgol_filter(B_R, B_window, B_order, mode='wrap')
B_Z = savgol_filter(B_Z, B_window, B_order, mode='wrap')
B_R = savgol_filter(B_R, B_window, B_order, mode='wrap')
B_Z = savgol_filter(B_Z, B_window, B_order, mode='wrap')
B_R[0] = B_R[-1]
B_Z[0] = B_Z[-1]
# plt.figure()
# plt.plot(B_R)
# plt.title("after")

B = np.sqrt((B_R**2.0) + (B_Z**2.0))
B_2 = (B_R**2.0) + (B_Z**2.0)




#checking to make sure that upsampling hasn't messed with periodicity. 
print("L dense periodicity check: " + str(L-l_dense[-1]))
print("R dense periodicity check: " + str((R_dense[-1] - R_dense[0])/(R_dense[0])))
print("Z dense periodicity check: " + str((Z_dense[-1] - Z_dense[0])/Z_dense[0]))


#Curvature Calculation _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- 

#Getting the derivatives that I need for curvature. 
dr_window = N_dense/10
dr_order = 3


spline_R_dense = CubicSpline(l_dense, R_dense, bc_type='periodic')
spline_Z_dense = CubicSpline(l_dense, Z_dense, bc_type='periodic')

R1_dense = spline_R_dense.derivative(1)(l_dense)
Z1_dense = spline_Z_dense.derivative(1)(l_dense)
R1_dense = savgol_filter(R1_dense, dr_window, dr_order, mode='wrap')
Z1_dense = savgol_filter(Z1_dense, dr_window, dr_order, mode='wrap')
R1_dense = savgol_filter(R1_dense, dr_window, dr_order, mode='wrap')
Z1_dense = savgol_filter(Z1_dense, dr_window, dr_order, mode='wrap')
R1_dense[-1] = R1_dense[0]
Z1_dense[-1] = Z1_dense[0]

spline_R1_dense = CubicSpline(l_dense, R1_dense, bc_type='periodic')
spline_Z1_dense = CubicSpline(l_dense, Z1_dense, bc_type='periodic')

R2_dense = spline_R1_dense.derivative(1)(l_dense)
Z2_dense = spline_Z1_dense.derivative(1)(l_dense)
# plt.figure()
# plt.plot(R2_dense)
# plt.title("before")
R2_dense = savgol_filter(R2_dense, dr_window, dr_order, mode='wrap')
Z2_dense = savgol_filter(Z2_dense, dr_window, dr_order, mode='wrap')
R2_dense = savgol_filter(R2_dense, dr_window, dr_order, mode='wrap')
Z2_dense = savgol_filter(Z2_dense, dr_window, dr_order, mode='wrap')

R2_dense[-1] = R2_dense[0]
Z2_dense[-1] = Z2_dense[0]

# plt.figure()
# plt.plot(R2_dense)
# plt.title("after")




#tangent vector is \vec{R1_dense, Z1_dense}
tangent_window = N_dense/10
tangent_order = 2

unit_tangent = np.sqrt((R1_dense**2.0)+(Z1_dense**2.0))
normal_mag = np.sqrt((R2_dense**2.0)+(Z2_dense**2.0))
unit_normal_R = R2_dense/normal_mag
unit_normal_Z = Z2_dense/normal_mag


# plt.figure()
# plt.plot(unit_tangent)
# plt.title("before")
unit_normal_R = savgol_filter(unit_normal_R, tangent_window, tangent_order, mode='wrap')
unit_normal_Z = savgol_filter(unit_normal_Z, tangent_window, tangent_order, mode='wrap')
unit_normal_R[-1] = unit_normal_R[0]
unit_normal_Z[-1] = unit_normal_Z[0]
# plt.figure()
# plt.plot(unit_tangent)
# plt.title("after")



#curvature. _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
kappa_dense = ((R1_dense*Z2_dense) - (Z1_dense*R2_dense)) / ((R1_dense**2) + (Z1_dense**2))**1.5
kappa_dense_R = kappa_dense*unit_normal_R
kappa_dense_Z = kappa_dense*unit_normal_Z
# plt.figure()
# plt.plot(kappa_dense)
# plt.title("before")

kappa_dense = savgol_filter(kappa_dense, window_length, polyorder, mode='wrap')
kappa_dense_R = savgol_filter(kappa_dense_R, window_length, polyorder, mode='wrap')
kappa_dense_Z = savgol_filter(kappa_dense_Z, window_length, polyorder, mode='wrap')

kappa_dense[-1] = kappa_dense[0]
kappa_dense_R[-1] = kappa_dense_R[0]
kappa_dense_Z[-1] = kappa_dense_Z[0]
# plt.figure()
# plt.plot(kappa_dense)
# plt.title("after")



#grad psi dot kappa _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
grad_psi_dot_kappa = (kappa_dense_R*dpsi_dR_dense) + (kappa_dense_Z*dpsi_dR_dense)


#d(grad_psi)/dl _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


d_grad_psi_dl = spline_grad_psi.derivative(1)(l_dense)

closest_index = 0
closest_distance = 100
B_min = 100000
B_test = 100000
#Finding the spot closest to the surface of the ring:
for i_l in range(0, len(l_dense)):
    #B_test = B[i_l]
    d = np.sqrt(((coil_major_radius - R_dense[i_l])**2.0) + ((coil_height - Z_dense[i_l])**2.0))
    # if B_test < B_min:
    #     B_min = B_test
    #     closest_index = i_l
    if d < closest_distance:
        closest_distance = d
        closest_index = i_l
        
        
print("Min field " + str(B_min))
print("at index "+str(closest_index))
print("(R, Z) = (" + str(R_dense[closest_index]) + ", " + str(Z_dense[closest_index]) +")")
#Rearranging all data so that this closest index is the start, so that I can set ICs there and integrate along. 

R_dense = swap_order(R_dense, closest_index)
R1_dense = swap_order(R1_dense, closest_index)
R2_dense = swap_order(R2_dense, closest_index)
Z_dense = swap_order(Z_dense, closest_index)
Z1_dense = swap_order(Z1_dense, closest_index)
Z2_dense = swap_order(Z2_dense, closest_index)
B = swap_order(B, closest_index)
B_2 = swap_order(B_2, closest_index)
grad_psi = swap_order(grad_psi, closest_index)
grad_psi_2 = swap_order(grad_psi_2, closest_index)
kappa_dense = swap_order(kappa_dense, closest_index)
kappa_dense_R = swap_order(kappa_dense_R, closest_index)
kappa_dense_Z = swap_order(kappa_dense_Z, closest_index)
grad_psi_dot_kappa = swap_order(grad_psi_dot_kappa, closest_index)
dpsi_dR_dense = swap_order(dpsi_dR_dense, closest_index)
dpsi_dZ_dense = swap_order(dpsi_dZ_dense, closest_index)


#Reinforcing Periodicity

R_dense[-1] = R_dense[0]
R1_dense[-1] = R1_dense[0]
R2_dense[-1] = R2_dense[0]
Z_dense[-1] = Z_dense[0]
Z1_dense[-1] = Z1_dense[0]
Z2_dense[-1] = Z2_dense[0]
B[-1] = B[0]
B_2[-1] = B_2[0]
grad_psi[-1] = grad_psi[0]
grad_psi_2[-1] = grad_psi_2[0]
kappa_dense[-1] = kappa_dense[0]
kappa_dense_R[-1] = kappa_dense_R[0]
kappa_dense_Z[-1] = kappa_dense_Z[0]
grad_psi_dot_kappa[-1] = grad_psi_dot_kappa[0]
dpsi_dR_dense[-1] = dpsi_dR_dense[0]
dpsi_dZ_dense[-1] = dpsi_dZ_dense[0]



# plt.figure()
# plt.plot(dpsi_dR_dense)
# plt.title("before")

grad_psi_R_spline = CubicSpline(l_dense, dpsi_dR_dense, bc_type='periodic')
grad_psi_Z_spline = CubicSpline(l_dense, dpsi_dR_dense, bc_type='periodic')
dgrad_psi_R_l = grad_psi_R_spline.derivative(1)(l_dense)
dgrad_psi_Z_l = grad_psi_Z_spline.derivative(1)(l_dense)
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
    return(2*((dpsi_dR_dense[l]*dgrad_psi_R_l[l])+(dpsi_dZ_dense[l]*dgrad_psi_Z_l[l])))

def A_3(l):
    return(((2*mu_0)/B_2[l])*(grad_psi_2[l])*(dp_dpsi)*(grad_psi_dot_kappa[l]))


coeff_window = N_dense/10
coeff_order = 2
A_1_array = np.zeros(N_dense)
A_2_array = np.zeros(N_dense)
A_3_array = np.zeros(N_dense)


for l in range(0, N_dense):
    A_1_array[l] = A_1(l)
    A_2_array[l] = A_2(l)
    A_3_array[l] = A_3(l)

    
A_1_array = savgol_filter(A_1_array, coeff_window, coeff_order, mode='wrap')    
A_2_array = savgol_filter(A_2_array, coeff_window, coeff_order, mode='wrap')    
A_3_array = savgol_filter(A_3_array, coeff_window, coeff_order, mode='wrap')    
#Enforcing periodicty by hand. 
A_1_array[0] = A_1_array[-1]
A_2_array[0] = A_2_array[-1]
A_3_array[0] = A_3_array[-1]
A1 = CubicSpline(l_dense, A_1_array, bc_type='periodic')
A2 = CubicSpline(l_dense, A_2_array, bc_type='periodic')
A3 = CubicSpline(l_dense, A_3_array, bc_type='periodic')


def rhs(l, y):
    y1, y2 = y
    return [y2, -(A2(l)*y2 + A3(l)*y1)/A1(l)]

solution = solve_ivp(rhs, t_span=(l_dense[0], l_dense[-1]), y0=[0.1, 0.0001], dense_output=True)
gamma = solution.sol(l_dense)
gamma = gamma[0]
print(len(gamma))
#Now we can convert this to the real space displacement of the lowest order pertubation. 
displacement = np.zeros(N_dense)
def disp_n(l, n_mode=10):
    return((n_mode*grad_psi[l]/B[l])*gamma[l])
    
for l in range(0, N_dense):
    displacement[l] = disp_n(l)
    
victim = np.zeros(N_dense)
for l in range(0, N_dense):
    victim[l] = (((2*mu_0)/B_2[l])*(grad_psi_2[l])*(dp_dpsi)*(grad_psi_dot_kappa[l]))
    
        
# #Hunting divergences. 
# plt.figure()
# plt.plot(l_dense, victim)
# plt.title("Victim")

#Plotting the contour. 
# plt.figure()
# plt.plot(R_dense, Z_dense)
# plt.plot(R_dense[0], Z_dense[0], 'rx')
# plt.title("Contour in Question")

plt.figure()
plt.plot(l_dense, displacement)
plt.xlabel("l")
plt.ylabel("pertubation size")
plt.title("First order mode shape")

#Component Checks
# plt.figure()
# plt.plot(l_dense, dpsi_dR_dense)
# plt.xlabel("l")
# plt.ylabel("dpsi_dR_dense")

# plt.figure()
# plt.plot(l_dense, dpsi_dZ_dense)
# plt.xlabel("l")
# plt.ylabel("dpsi_dZ_dense")




# plt.figure()
# plt.plot(l_dense, grad_psi_dot_kappa)
# plt.xlabel("l")
# plt.ylabel("grad psi dot kappa")

plt.figure()
plt.plot(l_dense, B_2)
plt.xlabel("l")
plt.ylabel("magnetic field squared")

# plt.figure()
# plt.plot(l_dense, grad_psi_2)
# plt.xlabel("l")
# plt.ylabel("grad psi field squared")


# #Coefficient Checks
plt.figure()
plt.plot(l_dense, A_1_array)
plt.xlabel("l")
plt.ylabel("A_1")
plt.title("Checking for zeroes of first coefficient")

plt.figure()
plt.plot(l_dense, A_2_array)
plt.xlabel("l")
plt.ylabel("A_2")
plt.title("Checking for zeroes of second coefficient")

plt.figure()
plt.plot(l_dense, A_3_array)
plt.xlabel("l")
plt.ylabel("A_3")
plt.title("Checking for zeroes of third coefficient")



plt.figure()
plt.plot(l_dense, gamma) 
plt.xlabel('l')
plt.ylabel('Gamma(l)')
plt.title('Ballooning Eigenvalue')
plt.grid(True)



plt.show()