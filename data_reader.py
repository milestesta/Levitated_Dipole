import numpy as np
import matplotlib.pyplot as plt

# Load the data from text files
pressure = np.loadtxt("pressure.txt")
psi = np.loadtxt("psi.txt")
label = np.loadtxt("labels.txt")
gs_check = np.loadtxt("GS_check.txt")
source_check = np.loadtxt("source_grid.txt")


# Get grid dimensions
NR, NZ = pressure.shape
print(NR)
print(NZ)
# Create grid coordinates
R = np.linspace(0, 1, NR)  
Z = np.linspace(0, 1, NZ)

# Create meshgrid
R_grid, Z_grid = np.meshgrid(R, Z, indexing='ij')


validation_ratio = (np.abs((gs_check-source_check)/source_check))

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

plt.figure(figsize=(8, 6))
plt.contourf(R_grid, Z_grid, label, levels = 6)#levels=np.linspace(0, 10, 50))
plt.colorbar(label="Grid Labels")
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Plotting the gridpoint labels")

plt.figure(figsize=(8, 6))
plt.contourf(R_grid, Z_grid, gs_check, levels = 6)#levels=np.linspace(0, 10, 50))
plt.colorbar(label="GS values")
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Plotting the solution check")

plt.figure(figsize=(8, 6))
plt.contourf(R_grid, Z_grid, source_check, levels = 6)#levels=np.linspace(0, 10, 50))
plt.colorbar(label="source values")
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Plotting the solution check")

plt.figure(figsize=(8, 6))
plt.contourf(R_grid, Z_grid, validation_ratio, levels = np.linspace(0, 1, 50))#levels=np.linspace(0, 10, 50))
plt.colorbar(label="ratio of source to solution")
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Plotting the ratio of source to solution")


plt.show()
