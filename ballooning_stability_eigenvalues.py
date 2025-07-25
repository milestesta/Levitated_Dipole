import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, RectBivariateSpline, CubicSpline, UnivariateSpline
from scipy.signal import savgol_filter
from scipy.constants import mu_0
from scipy.integrate import solve_bvp, solve_ivp
from skimage.measure import find_contours


N_levels = 10

#coil parameters
coil_height = 1.5
coil_major_radius = 0.4
coil_minor_radius = 0.1


pressure_max = 689476.0
psi_coil = 0.000501772
psi_max = 2.50886

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





# # Plot psi contour
# plt.figure(figsize=(8, 6))
# CS = plt.contour(R_dense, Z_dense, psi_dense, levels = 50)#np.linspace(0, 0.5, 50))
# plt.colorbar(label="Psi Dense")
# # plt.clabel(CS, inline=1, fontsize=10)
# plt.xlabel("R")
# plt.ylabel("Z")
# plt.title("Contour Plot of Psi")
# plt.show()


#function definitions. _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
def pressure_value(flux):
    return((pressure_max*(1 - np.cos((flux/psi_max)*np.pi))))
    
def pressure_derivative(flux):
    pressure_max = 1e6
    if flux < psi_coil:
        return(0)
    elif flux >= psi_coil:
        psi_normal = ((flux - psi_coil)/(psi_max - psi_coil))
        conversion = psi_max-psi_coil
        return(pressure_max*conversion*np.pi*(np.sin(np.pi*psi_normal)))
    else:
        return(0)

def swap_order(X, I):
    #X is the array we want to swap the order of
    #I is the index of the element that we want to be first. 
    A = np.asarray(X[I:], float)
    B = np.asarray(X[0:I], float)
    Y = np.concatenate((A, B))
    return(Y)

def s_finder(y):
    #For estimating the filter level that I chuck into univariate. 
    dy = np.diff(y)
    s = ((np.std(dy)/np.sqrt(2))**2.0)*len(y)
    return(s)

def is_closed(set_of_points):
    #returns false for open, true for closed. 
    N_set = len(set_of_points)
    tol = ((DR_dense**2.0) + (DZ_dense**2.0))
    for n in range(0, N_set-1):
        R1 = set_of_points[n, 0]
        R2 = set_of_points[n+1, 0]
        Z1 = set_of_points[n, 1]
        Z2 = set_of_points[n+1, 1]
        D = ((R2-R1)**2.0) + ((Z2-Z1)**2.0)
        if D >= tol:
            return(False)
    print("aha!")
    return(True)


def array_zero_finder(domain, A):
    #given an array A, find where it crosses zero. 
    N_A = len(A)
    roots = []
    root_indices = []
    for i in range(1, N_A):
        if A[i]*A[i-1] < 0:
            roots.append((domain[i] + domain[i-1])/2.0)
            root_indices.append(int(i))
    return(roots, root_indices)

#Picking out the closed contours: _-____--___-__-__-




d_psi = (psi_dense.max())/N_levels


psi_levels = []
contours = []



def find_closed_contours_skimage(psi_grid, level, R_grid, Z_grid):
    contours = find_contours(psi_grid, level)
    closed = []
    levels = []
    for c in contours:
        if np.allclose(c[0], c[-1], atol=1e-8):
            # Convert from indices to R,Z
            r_coords = np.array(c[:, 1]*DR_dense)
            z_coords = np.array(c[:, 0]*DZ_dense)
            
            closed.append(np.column_stack([r_coords, z_coords]))
            levels.append(level)
    return closed, level



for n in range(0, N_levels):
    psi_level = (d_psi*n)
    cont, lev = find_closed_contours_skimage(psi_dense, psi_level, R_dense, Z_dense)
    n_at_at_level = len(cont)
    for k in range(0, n_at_at_level):
        psi_levels.append(lev)
        contours.append(cont[k])

N_closed = len(psi_levels)
print("found " + str(N_closed) + " out of " + str(N_levels) + " to be closed.")




def contour_eigenvalue_analysis(contour_points, contour_psi):
    R_contour, Z_contour = contour_points[:, 0], contour_points[:, 1] #R, Z values of each point on the contour. 
    R_indices = (R_contour/DR_dense) #Needed to get the field points along the contour. 
    Z_indices = (Z_contour/DZ_dense)    
    N = len(R_indices) #Number of points on the contour. 
    dp_dpsi = pressure_derivative(contour_psi)
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
    if grad_psi_window <  grad_psi_poly_order:
        grad_psi_window = grad_psi_poly_order+1


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
    
    if dgrad_dl_window < d_grad_dl_poly_order:
        dgrad_dl_window = d_grad_dl_poly_order+1
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

    # A1_arr = np.zeros(N)
    # A2_arr = np.zeros(N)
    # A3_arr = np.zeros(N)

    # for i in range(0, N): 
    #     A1_arr[i] = B[i]/mag_grad_psi_2[i]
    #     A2_arr[i] = 2.0*(d_psi_d_R_contour[i]*d_grad_psi_dl_R[i] + d_psi_d_Z_contour[i]*d_grad_psi_dl_Z[i])
    #     A3_arr[i] = (2.0*mu_0 / B_2[i]) * mag_grad_psi_2[i] * dp_dpsi * grad_psi_dot_kappa[i]




    # # plt.figure()
    # # plt.plot(l, A1_arr)
    # # plt.title("A1")
    # # plt.figure()
    # # plt.plot(l, A2_arr)
    # # plt.title("A2")
    # # plt.figure()
    # # plt.plot(l, A3_arr)
    # # plt.title("A3")

    # A1 = UnivariateSpline(l, A1_arr, s=s_finder(A1_arr))
    # A2 = UnivariateSpline(l, A2_arr, s=s_finder(A2_arr))
    # A3 = UnivariateSpline(l, A3_arr, s=s_finder(A3_arr))

        #Ode coefficients calculation: _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
    A1_array = B/mag_grad_psi_2
    A1 = UnivariateSpline(l, A1_array, s=s_finder(A1_array))
    A2_array = A1.derivative()(l)
    A2 = UnivariateSpline(l, A2_array, s=s_finder(A2_array))
    A3_array = (mu_0*dp_dpsi)*(1/(mag_grad_psi_2*B))*(grad_psi_dot_kappa)
    A3 = UnivariateSpline(l, A3_array, s=s_finder(A3_array))
    
    # plt.figure()
    # plt.plot(l, A1(l))
    # plt.title("A1")
    # plt.figure()
    # plt.plot(l, A2(l))
    # plt.title("A2")
    # plt.figure()
    # plt.plot(l, A3(l))
    # plt.title("A3")

    # A1_p = A1.derivative()(l)
    # A2_p = A2.derivative()(l)


    y_guess = np.zeros((2, N))
    # small sinusoidal perturbation
    y_guess[0] = 1e-7*np.cos(6*np.pi*l / l[-1])
    y_guess[1] = -(6*np.pi/l[-1])*1e-7*np.sin(6*np.pi*l / l[-1])

    # y_guess[0] = np.ones(N)*1e-5
    # y_guess[1] = np.zeros(N)


    # def rhs(l, y):
    #     # y[0] = Gamma, y[1] = dGamma/dl I keep messing this up
    #     return np.vstack([y[1], -(A2(l)*y[1] + A3(l)*y[0]) / A1(l)])
    
    def fun(x, Y):
        y, yp, lam = Y                # each of shape (M,)
        # A   = A_spline(x)             # → (M,)
        # dA  = dA_spline(x)            # → (M,)
        # B   = B_spline(x)             # → (M,)
        # W   = W_spline(x) if W_spline else 1.0  # → scalar or (M,)

        dy0 = yp
        dy1 = - (A2(x)/A1(x))*yp - (A3(x) + lam)/A1(x) * y

        # THIS is the key:
        dp  = np.zeros_like(x)        # → (M,)

        return np.vstack([dy0, dy1, dp])   # now shape (3, M)
    # def rhs(l, y):
    #     y1, y2 = y
    #     return [y2, -(A2(l)*y2 + A3(l)*y1)/A1(l)]

    # 2) periodic BCs: y(0)=y(L), y'(0)=y'(L) REMEMBER WHICH INDEX CORRESPONDS TO DERIVATIVE YOU MORON. 
    # def bc(ya, yb):
    #     return np.array([ya[0] - yb[0], ya[1] - yb[1]])
    
    def bc(Y0, YL):
    # Y0 = [y0(0), y1(0), p]
    # YL = [y0(L), y1(L), p]
        return np.array([
            Y0[0] - YL[0],   # y(0) - y(L) = 0
            Y0[1] - YL[1],   # y'(0) - y'(L) = 0
            Y0[0] - 1.0      # y(0) = 1 normalization
        ])

    # x = np.linspace(0, L, 50)

    # initial guess: y = cos(2πx/L), y' = -2π/L sin(2πx/L), p = 0.1
    Y_init = np.vstack([
        np.cos(2*np.pi*l/L),
    -2*np.pi/L * np.sin(2*np.pi*l/L),
        0.1 * np.ones_like(l)
    ])
    
    sol = solve_bvp(fun, bc, l, Y_init, tol=1e-2, max_nodes=10*N)

    if sol.status != 0:
        raise RuntimeError("BVP solver did not converge")

    # eigenvalue and eigenfunction:
    lambda_n = sol.y[2,0]
    gamma      = sol.y[0]       # on sol.x mesh

    print("Found eigenvalue λ = " + str(lambda_n) + "")
    # ivp_sol = solve_ivp(rhs, t_span=(0, L), y0=[1e-7, 1e-3], dense_output=True) 
    # solution = solve_bvp(rhs, bc, l, y_guess, tol = 1e-10, max_nodes = 100)
    # [gamma, dgamma] = ivp_sol.sol(l)
    # plt.plot(sol.x, gamma, '-o')
    # plt.xlabel("x")
    # plt.ylabel("y(x)")
    # plt.title(f"Eigenmode, λ = {lambda_n:.4f}")
    # plt.show()
    displacement = np.zeros(N)
    def disp_n(i, n_mode=10):
        return((n_mode*mag_grad_psi[i]/B[i])*gamma[i])
    l_solution = np.linspace(0, L, len(gamma))        
    for i in range(0, N):
        displacement[i] = disp_n(i)
    root, root_ind = array_zero_finder(l_solution, gamma)

    return(l_solution, gamma, displacement, root, root_ind, lambda_n)

contour_lengths = []
gammas = []
displacements = []
roots = []
root_indices = []
is_unstable = np.zeros(N_closed)
eigenfunctions = []


for n in range(0, N_closed):
    if (len(contours[n]) > 100):
        ll, gam, disp, rot, rot_ind, eigen = contour_eigenvalue_analysis(contours[n], psi_levels[n])
        contour_lengths.append(ll)
        gammas.append(gam)
        displacements.append(disp)
        roots.append(rot)
        root_indices.append(rot_ind)
        eigenfunctions.append(eigen)
        if (len(rot) != 0):
            is_unstable[n] = 1

if np.sum(is_unstable) > 0:
    print("The Equilibrium was Ballooning Unstable!")
else:
    print("The Equilibrium was Ballooning Stable!")
    
print("_-_-_-_-_-_-_-_-_-_-")
print("Starting the prints and plots...")
    
for n in range(0, N_closed):
    #print(len(contours[n]))
    if is_unstable[n] == 1:
        contour = contours[n]
        R_set = contour[:, 0]
        Z_set = contour[:, 1]
        plt.plot(R_set, Z_set, "r")
    if is_unstable[n] == 0:
        contour = contours[n]
        R_set = contour[:, 0]
        Z_set = contour[:, 1]
        plt.plot(R_set, Z_set, "g")      
plt.title("All Contours")
plt.figure()

for n in range(0, N_closed):
    if is_unstable[n] == 1:
        contour = contours[n]
        R_set = contour[:, 0]
        Z_set = contour[:, 1]
        plt.plot(R_set, Z_set, "r")
plt.title("Unstable Contours")
plt.figure()

for n in range(0, N_closed):
    if is_unstable[n] == 0:
        contour = contours[n]
        R_set = contour[:, 0]
        Z_set = contour[:, 1]
        plt.plot(R_set, Z_set, "g")
plt.title("Stable Contours")
print(is_unstable)

plt.figure()
plt.title("A representative 'solution'")
n_roots = len(roots[1])
rot_plot = np.zeros(n_roots)
plt.plot(contour_lengths[1], gammas[1])
plt.plot(roots[1], rot_plot, "rx")
     
# R_contour, Z_contour = contour_points[:, 0], contour_points[:, 1] #R, Z values of each point on the contour. 
# R_indices = (R_contour/DR_dense) #Needed to get the field points along the contour. 
# Z_indices = (Z_contour/DZ_dense)
# N = len(R_indices) #Number of points on the contour. 
# print("we have " + str(N) + " points on the contour.")

# pressure_max = 689476.0
# psi_max = 1.0
print(eigenfunctions)
plt.show()