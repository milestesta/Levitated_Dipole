import numpy as np
import matplotlib.pyplot as plt

# Load the data from text files
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

print(NR)
print(NZ)
print(Z)

R_grid, Z_grid = np.meshgrid(R, Z, indexing='ij')


validation_ratio = (np.abs(gs_check - source_check)/np.abs(gs_check))*100 #percentage error. 
# Plot psi contour
plt.figure(figsize=(8, 6))
CS = plt.contour(R_grid, Z_grid, psi, levels=50)#np.linspace(0, 0.5, 50))
plt.colorbar(label="Psi")
# plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Contour Plot of Psi")


# Plot pressure contour
plt.figure(figsize=(8, 6))
plt.contour(R_grid, Z_grid, pressure, levels = 50)#levels=np.linspace(0, 10, 50))
plt.colorbar(label="Pressure")
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Contour Plot of Pressure")

# plt.figure(figsize=(8, 6))
# plt.contour(R_grid, Z_grid, label, levels = 6)#levels=np.linspace(0, 10, 50))
# plt.colorbar(label="Grid Labels")
# plt.xlabel("R")
# plt.ylabel("Z")
# plt.title("Plotting the gridpoint labels")

plt.figure(figsize=(8, 6))
plt.contour(R_grid, Z_grid, gs_check, levels=50)
plt.colorbar(label="GS Operator")
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Plotting the solution check")

plt.figure(figsize=(8, 6))
plt.contour(R_grid, Z_grid, source_check, levels=50)#levels=np.linspace(0, 10, 50))
plt.colorbar(label="source values")
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Plotting the solution check")

# plt.figure(figsize=(8, 6))
# plt.contour(R_grid, Z_grid, validation_ratio, levels = 50)#levels=np.linspace(0, 10, 50))
# plt.colorbar(label="difference between source and solution")
# plt.xlabel("R")
# plt.ylabel("Z")
# plt.title("source vs GS acting on psi")

plt.figure(figsize=(8, 6))
plt.contour(R_grid, Z_grid, dipole_field, levels = 50)#levels=np.linspace(0, 10, 50))
plt.colorbar(label="Pre-plasma field")
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Plot of psi before we add the plasma")

plt.figure(figsize=(8, 6))
plt.contourf(R_grid, Z_grid, validation_ratio, levels = np.linspace(0, 10, 50))#levels=np.linspace(0, 10, 50))
plt.colorbar(label="Error")
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Percentage error in the solution")


plt.show()
