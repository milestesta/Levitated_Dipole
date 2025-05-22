import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, RectBivariateSpline, CubicSpline, UnivariateSpline
from scipy.signal import savgol_filter
from scipy.constants import mu_0
from scipy.integrate import solve_bvp, solve_ivp

#coil parameters
coil_height = 1.5
coil_major_radius = 0.4
coil_minor_radius = 0.1




# Loading the data from text files_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
pressure = np.loadtxt("pressure.txt")
psi = np.loadtxt("psi.txt")
label = np.loadtxt("labels.txt")
gs_check = np.loadtxt("GS_check.txt")
source_check = np.loadtxt("source_grid.txt")
dipole_field = np.loadtxt("initial_grid.txt")
R = np.loadtxt("R_grid.txt")
Z = np.loadtxt("Z_grid.txt")
# Get grid dimensions
NR = len(R)
NZ = len(Z)
n_R_upsample = 1000
n_Z_upsample = 1000
# print(NR)
# print(NZ)
# print(Z)

#Upsampling the overall Psi grid. I feel like one big round of upsampling is better than lots of little ones. -_-_-_-_-_-_-_-_
R_grid, Z_grid = np.meshgrid(R, Z, indexing='ij')
R_min = np.min(R_grid)
Z_min = np.min(Z_grid)
Z_max = np.max(Z_grid)
R_max = np.max(R_grid)

R_sparse = np.linspace(R_min, R_max, NR)
Z_sparse = np.linspace(Z_min, Z_max, NZ)

psi_dense_interp = RectBivariateSpline(R_sparse, Z_sparse, psi.T)

R_dense = np.linspace(R[0], R[-1], n_R_upsample)
Z_dense = np.linspace(Z[0], Z[-1], n_Z_upsample)

psi_dense = psi_dense_interp(R_dense, Z_dense)

DR_dense  = (R[-1] - R[0])/n_R_upsample
DZ_dense  = (Z[-1] - Z[0])/n_Z_upsample


psi_level = (psi_dense.max() - psi_dense.min())*0.61#psi value of the contour I care about. 
fig, ax = plt.subplots()
contour_1 = ax.contour(R_dense, Z_dense, psi_dense, levels=[psi_level])
contour_path = contour_1.collections[0].get_paths()
contour_points = contour_path[0].vertices  # (N, 2) array of (R, Z) points
R_contour, Z_contour = contour_points[:, 0], contour_points[:, 1] #R, Z values of each point on the contour. 
R_indices = (R_contour/DR_dense) #Needed to get the field points along the contour. 
Z_indices = (Z_contour/DZ_dense)
N = len(R_indices) #Number of points on the contour. 
print("we have " + str(N) + " points on the contour.")

pressure_max = 689476.0
psi_max = 1.0



#function definitions. _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
def pressure_value(flux):
    return((pressure_max*(1 - np.cos((flux/psi_max)*np.pi))))
    
def pressure_derivative(flux):
    return((pressure_max*np.pi/psi_max)*(np.sin((flux)*np.pi/psi_max)))

def swap_order(X, I):
    #X is the array we want to swap the order of
    #I is the index of the element that we want to be first. 
    
    N = len(X)
    A = np.asarray(X[I:], float)
    B = np.asarray(X[0:I], float)
    Y = np.concatenate((A, B))
    return(Y)

def s_finder(y):
    #For estimating the filter level that I chuck into univariate. 
    dy = np.diff(y)
    s = ((np.std(dy)/np.sqrt(2))**2.0)*len(y)
    return(s)


dp_dpsi = pressure_derivative(psi_level)
print("pressure gradient: " + str(dp_dpsi))

#Calculating arc length along the contour. -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
l = np.zeros(N)
L = 0
for n in range(1, N):
    L += np.sqrt(((R_contour[n] - R_contour[n-1])**2.0) + ((Z_contour[n] - Z_contour[n-1])**2.0))
    l[n] = L

print("Contour has length L = " + str(L))


#Calculating the derivatives of R and Z on the contour. Needed for curvature. -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

#Spline fitting to make this neater. 
R_spline = CubicSpline(l, R_contour, bc_type='periodic')
Z_spline = CubicSpline(l, Z_contour, bc_type='periodic')

R_univ = UnivariateSpline(l, R_contour, s=1e-5)
# plt.figure()
# plt.title("Comparing cubic an univariate")
# plt.plot(l, R_spline(l), "r")
# plt.figure()
# plt.plot(l, R_univ(l), "b")

R_points_splined = R_spline(l)
Z_points_splined = Z_spline(l)

dR_1 = R_spline.derivative(1)(l)
dZ_1 = Z_spline.derivative(1)(l)

dR_1[0] = dR_1[-1]
dZ_1[0] = dZ_1[-1]

#Have to smooth the first derivatives again or I basically just get noise.
#Sav-Gol wasn't working, so did some digging an univ seemed to work better. 
# plt.figure()
# plt.plot(l, dZ_1)
# plt.title("Dr1 before")

dR_spline = UnivariateSpline(l, dR_1, s=s_finder(dR_1))
dZ_spline = UnivariateSpline(l, dZ_1, s=s_finder(dZ_1))

# plt.figure()
# plt.plot(l, dZ_spline(l))
# plt.title("Dr1 after")
dR_2 = dR_spline.derivative()(l)
dZ_2 = dZ_spline.derivative()(l)


dR_2[0] = dR_2[-1]
dZ_2[0] = dZ_2[-1]

# plt.figure()
# plt.plot(l, dR_2)
# plt.title("Dr2 before")
#Now I have to make everything periodic again....

# derivative_smoothing_window = N//30 + (1 - (N//30)%2) 
# derivative_smoothing_poly_order = 3

# dR_1 = savgol_filter(dR_1, derivative_smoothing_window, derivative_smoothing_poly_order, mode='wrap')
# dZ_1 = savgol_filter(dZ_1, derivative_smoothing_window, derivative_smoothing_poly_order, mode='wrap')
# dR_2 = savgol_filter(dR_2, derivative_smoothing_window, derivative_smoothing_poly_order, mode='wrap')
# dZ_2 = savgol_filter(dZ_2, derivative_smoothing_window, derivative_smoothing_poly_order, mode='wrap')




#Calculating curvature:-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

unit_tangent = np.sqrt((dR_1**2.0)+(dZ_1**2.0))

normal_mag = np.sqrt((dR_2**2.0)+(dZ_2**2.0))
unit_normal_R = dR_2/normal_mag
unit_normal_Z = dZ_2/normal_mag

kappa = ((dR_1*dZ_2) - (dZ_1*dR_2)) / ((dR_1**2) + (dZ_1**2))**1.5
kappa_R = kappa*unit_normal_R
kappa_Z = kappa*unit_normal_Z

# plt.figure()
# plt.title("Curvature Magnitude")
# plt.plot(l, kappa_R)




#Calculating Grad_psi on the contour. -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
d_psi_d_Z, d_psi_d_R = np.gradient(psi_dense, DZ_dense, DR_dense)
d_psi_d_R_contour = np.zeros(N) 
d_psi_d_Z_contour = np.zeros(N)
for point in range(0, N):
    d_psi_d_R_contour[point] = d_psi_d_R[int(R_indices[point])][int(Z_indices[point])]
    d_psi_d_Z_contour[point] = d_psi_d_Z[int(R_indices[point])][int(Z_indices[point])]

#smoothing grad psi. 
grad_psi_window = N//10 + (1 - (N//10)%2) 
grad_psi_poly_order = 3


# plt.figure()
# plt.plot(l, d_psi_d_R_contour)
# plt.title('checking dpsidR before')



# d_psi_d_R_contour = savgol_filter(d_psi_d_R_contour, grad_psi_window, grad_psi_poly_order, mode='wrap')
# d_psi_d_Z_contour = savgol_filter(d_psi_d_Z_contour, grad_psi_window, grad_psi_poly_order, mode='wrap')
# d_psi_d_R_contour[0] = d_psi_d_R_contour[-1]
# d_psi_d_Z_contour[0] = d_psi_d_Z_contour[-1]


# plt.figure()
# plt.plot(l, d_psi_d_R_contour)
# plt.title('checking dpsidR after')
#I now have to calculate the derivatives of the components of grad psi w.r.t l. 

d_psi_d_R_spline = UnivariateSpline(l, d_psi_d_R_contour, s=s_finder(d_psi_d_R_contour))
d_psi_d_Z_spline = UnivariateSpline(l, d_psi_d_Z_contour, s=s_finder(d_psi_d_Z_contour))


d_grad_psi_dl_R = d_psi_d_R_spline.derivative()(l)
d_grad_psi_dl_Z = d_psi_d_Z_spline.derivative()(l)

# plt.figure()
# plt.plot(l, d_grad_psi_dl_R)
# plt.title('checking grad deriv. BEFORE')

#and making it periodic again. 

d_grad_psi_dl_R  = savgol_filter(d_grad_psi_dl_R, grad_psi_window, grad_psi_poly_order, mode='wrap')
d_grad_psi_dl_Z  = savgol_filter(d_grad_psi_dl_Z , grad_psi_window, grad_psi_poly_order, mode='wrap')
d_grad_psi_dl_R[0] = d_grad_psi_dl_R[-1]
d_grad_psi_dl_Z[0] = d_grad_psi_dl_Z[-1]

# plt.figure()
# plt.plot(l, d_grad_psi_dl_R)
# plt.title('checking grad deriv. AFTER')
mag_grad_psi = np.sqrt((d_psi_d_R_contour*d_psi_d_R_contour)+(d_psi_d_Z_contour*d_psi_d_Z_contour))
mag_grad_psi_2 = mag_grad_psi*mag_grad_psi



# plt.figure()
# plt.plot(l, mag_grad_psi)
# plt.title('|\grad psi|')

#grad psi dot kappa _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
grad_psi_dot_kappa = (kappa_R*d_psi_d_R_contour) + (kappa_Z*d_psi_d_Z_contour)




#d(grad_psi)/dl _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
grad_psi_univ = UnivariateSpline(l, mag_grad_psi, s=s_finder(mag_grad_psi))
new_mag_grad_psi = grad_psi_univ(l)
# plt.figure()
# plt.plot(l, new_mag_grad_psi)
# plt.title('new |\grad psi|')
d_grad_psi_dl = grad_psi_univ.derivative()(l)
d_grad_psi_dl[0] = d_grad_psi_dl[-1]


#smoothing: 
dgrad_dl_window = N//10 + (1 - (N//10)%2) 
d_grad_dl_poly_order = 3
d_grad_psi_dl = savgol_filter(d_grad_psi_dl, dgrad_dl_window, d_grad_dl_poly_order, mode='wrap')
# plt.figure()
# plt.title('dgrad_psi_dl')
# plt.plot(d_grad_psi_dl)

#Magnetic Field Calculation _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
B_R = (-1.0)*d_psi_d_Z_contour/R_contour
B_Z = d_psi_d_R_contour/R_contour



B = np.sqrt((B_R**2.0) + (B_Z**2.0))
B_2 = (B_R**2.0) + (B_Z**2.0)

# plt.figure()
# plt.plot(l, B_2)
# plt.title("B2")

#Swapping the order of the arrays to all have the start of the eigenmode as the start of the array. 
#For the moment, I assume that the eigenmode starts with zero magnitude (and finite gradient) at the point
#with the smallest magnetic field. 


closest_index = 0
closest_distance = 100
B_min = 100000
B_test = 100000
#Finding the spot closest to the surface of the ring:
for i_l in range(0, len(l)):
    B_test = B[i_l]
    # d = np.sqrt(((coil_major_radius - R_contour[i_l])**2.0) + ((coil_height - Z_contour[i_l])**2.0))
    if B_test < B_min:
        B_min = B_test
        closest_index = i_l
    # if d < closest_distance:
    #     closest_distance = d
    #     closest_index = i_l
        
        
print("Min field " + str(B_min))
print("at index "+str(closest_index))
print("(R, Z) = (" + str(R_contour[closest_index]) + ", " + str(Z_contour[closest_index]) +")")
#Rearranging all data so that this closest index is the start, so that I can set ICs there and integrate along. 

R_contour = swap_order(R_contour, closest_index)
d_R1 = swap_order(dR_1, closest_index)
d_R2 = swap_order(dR_2, closest_index)
Z_contour = swap_order(Z_contour, closest_index)
dZ_1 = swap_order(dZ_1, closest_index)
dZ_2 = swap_order(dZ_2, closest_index)
B = swap_order(B, closest_index)
B_2 = swap_order(B_2, closest_index)
mag_grad_psigrad_psi = swap_order(mag_grad_psi, closest_index)
mag_grad_psi_2 = swap_order(mag_grad_psi_2, closest_index)
kappa = swap_order(kappa, closest_index)
kappa_R = swap_order(kappa_R, closest_index)
kappa_Z = swap_order(kappa_Z, closest_index)
grad_psi_dot_kappa = swap_order(grad_psi_dot_kappa, closest_index)
d_psi_d_R_contour = swap_order(d_psi_d_R_contour, closest_index)
d_psi_d_Z_contour = swap_order(d_psi_d_Z_contour, closest_index)
d_grad_psi_dl_R = swap_order(d_grad_psi_dl_R, closest_index)
d_grad_psi_dl_Z = swap_order(d_grad_psi_dl_Z, closest_index)


#Reinforcing Periodicity

R_contour[-1] = R_contour[0]
dR_1[-1] = dR_1[0]
dR_2[-1] = dR_2[0]
Z_contour[-1] = Z_contour[0]
dZ_1[-1] = dZ_1[0]
dZ_2[-1] = dZ_2[0]
B[-1] = B[0]
B_2[-1] = B_2[0]
mag_grad_psi[-1] = mag_grad_psi[0]
mag_grad_psi_2[-1] = mag_grad_psi_2[0]
kappa[-1] = kappa[0]
kappa_R[-1] = kappa_R[0]
kappa_Z[-1] = kappa_Z[0]
grad_psi_dot_kappa[-1] = grad_psi_dot_kappa[0]
d_psi_d_R_contour[-1] = d_psi_d_R_contour[0]
d_psi_d_Z_contour[-1] = d_psi_d_Z_contour[0]
d_grad_psi_dl_R[-1] = d_grad_psi_dl_R[0]
d_grad_psi_dl_Z[-1] = d_grad_psi_dl_Z[0]


#We now have everything we need for the ballooning equation. 


A1_arr = np.zeros(N)
A2_arr = np.zeros(N)
A3_arr = np.zeros(N)

for i in range(0, N): 
    A1_arr[i] = mag_grad_psi_2[i]
    A2_arr[i] = 2.0*(d_psi_d_R_contour[i]*d_grad_psi_dl_R[i] + d_psi_d_Z_contour[i]*d_grad_psi_dl_Z[i])
    A3_arr[i] = (2.0*mu_0 / B_2[i]) * mag_grad_psi_2[i] * dp_dpsi * grad_psi_dot_kappa[i]


# plt.figure()
# plt.plot(l, A1_arr)
# plt.title("A1")
# plt.figure()
# plt.plot(l, A2_arr)
# plt.title("A2")
# plt.figure()
# plt.plot(l, A3_arr)
# plt.title("A3")

A1 = UnivariateSpline(l, A1_arr, s=s_finder(A1_arr))
A2 = UnivariateSpline(l, A2_arr, s=s_finder(A2_arr))
A3 = UnivariateSpline(l, A3_arr, s=s_finder(A3_arr))


plt.figure()
plt.plot(l, A1(l))
plt.title("A1")
plt.figure()
plt.plot(l, A2(l))
plt.title("A2")
plt.figure()
plt.plot(l, A3(l))
plt.title("A3")

# A1_p = A1.derivative()(l)
# A2_p = A2.derivative()(l)


y_guess = np.zeros((2, N))
# small sinusoidal perturbation
y_guess[0] = 1e-7*np.cos(6*np.pi*l / l[-1])
y_guess[1] = -(6*np.pi/l[-1])*1e-7*np.sin(6*np.pi*l / l[-1])

# y_guess[0] = np.ones(N)*1e-5
# y_guess[1] = np.zeros(N)


def rhs(l, y):
    # y[0] = Gamma, y[1] = dGamma/dl I keep messing this up
    return np.vstack([y[1], -(A2(l)*y[1] + A3(l)*y[0]) / A1(l)])

# 2) periodic BCs: y(0)=y(L), y'(0)=y'(L) REMEMBER WHICH INDEX CORRESPONDS TO DERIVATIVE YOU MORON. 
def bc(ya, yb):
    return np.array([ya[0] - yb[0], ya[1] - yb[1]])

solution = solve_bvp(rhs, bc, l, y_guess, tol = 1e-10, max_nodes = 100)
gamma = solution.sol(l)
gamma = gamma[0]
print(len(gamma))
plt.figure()
plt.title("Gamma")
plt.plot(l, gamma)


#Now we can convert this to the real space displacement of the lowest order pertubation. 
displacement = np.zeros(N)
def disp_n(i, n_mode=10):
    return((n_mode*mag_grad_psi[i]/B[i])*gamma[i])
    
for i in range(0, N):
    displacement[i] = disp_n(i)
plt.figure()
plt.title('Real Space Mode shape')
plt.plot(l, displacement)

# plt.figure(figsize=(8, 6))
# CS = plt.contour(R_dense, Z_dense, psi_dense, levels=100)#np.linspace(0, 0.5, 50))
# plt.colorbar(label="Psi")
# # plt.clabel(CS, inline=1, fontsize=10)
# plt.xlabel("R")
# plt.ylabel("Z")
# plt.title("Contour Plot of Psi")
# plt.figure()
# plt.imshow(psi_dense)



plt.show()
