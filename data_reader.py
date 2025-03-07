import numpy as np
import matplotlib.pyplot as plt

# Load the data from text files
pressure = np.loadtxt("pressure.txt")
psi = np.loadtxt("psi.txt")

# Get grid dimensions
NR, NZ = pressure.shape
print(NR)
print(NZ)
# Create grid coordinates
R = np.linspace(0, 1, NR)  
Z = np.linspace(0, 1, NZ)

# Create meshgrid
R_grid, Z_grid = np.meshgrid(R, Z, indexing='ij')

# Plot psi contour
plt.figure(figsize=(8, 6))
plt.contour(R_grid, Z_grid, psi, levels=50)#np.linspace(0, 0.5, 50))
plt.colorbar(label="Psi")
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
plt.show()
