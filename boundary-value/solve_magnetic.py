from classes import Lattice
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

SIDE_LEN = int(input("Choose the system's side length: "))
METHOD = input("Choose the evolution method (Jacobi/Gauss-Seidel/SOR): ").lower()
if METHOD == "sor":
	omega = 1.78 # fastest
else:
	omega = None
TOLERANCE = float(input("Choose the accuracy of the solution: "))
# choose current as one likes
current = np.zeros(shape=(SIDE_LEN, SIDE_LEN, SIDE_LEN))
MID = SIDE_LEN // 2
current[MID, MID, 1:-1] = 1

state = Lattice(SIDE_LEN, current)

convergence = False
while not convergence:
	diff = state.evolve(method=METHOD, omega=omega) # read method's docstring for deeper understanding
	if diff <= TOLERANCE:
		convergence = True

B_field_squared = state.B_field**2 # square E_field components
B_field_strength = np.sqrt(np.sum(B_field_squared, axis=3)) # takes magnitude of E_field
distance = np.arange(1, MID - 1) # distance from cell next to CENTER up to boundary (excluding it)
# lines depndent only on distance to charge
Az_line = state.lattice[MID + 1:-1, MID, MID]
B_field_line = B_field_strength[MID + 1:-1, MID, MID]

# get data at midplane cut
Az_midplane = state.lattice[:, :, MID]
Bx_midplane = state.B_field[:, :, MID, 0]
By_midplane = state.B_field[:, :, MID, 1]

# save data
os.makedirs("./datafiles", exist_ok=True)
np.savetxt("./datafiles/vector_potential.dat", Az_midplane, header="Vector Potential Lattice")
np.savetxt("./datafiles/B_field.dat", state.B_field[:, :, MID].reshape(-1, 2), header="Magnetic Field Lattice")

x = np.arange(SIDE_LEN)
y = np.arange(SIDE_LEN)
X, Y = np.meshgrid(x, y)
fig, ax = plt.subplots()
contour = ax.contourf(X, Y, Az_midplane, levels=50, cmap="RdBu_r") # creates the contour with 50 discretized levels
colbar = plt.colorbar(contour, ax=ax, fraction=.046, pad=.04)
colbar.set_label(r"Vector Potential ($A_z$)")
jump = 4 # reduces the number of arrows that will be seen in the E_field plot
# normalize the E_field to get appropriate arrow lengths
B_mag_midplane = np.hypot(Bx_midplane, By_midplane)
B_mag_midplane[B_mag_midplane == 0] = 1
Bx_norm = Bx_midplane / B_mag_midplane
By_norm = By_midplane / B_mag_midplane
ax.quiver(X[::jump, ::jump], Y[::jump, ::jump],
        Bx_norm[::jump, ::jump], By_norm[::jump, ::jump],
        color="white",
        pivot="mid",
        alpha=.8) # creates a 2D slice of the 3D E_field
ax.set_title("Cut through midplane: Vector Potential and Magnetic Field")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal") # ensures the image doesn't get stretched
plt.savefig(f"./plots/magnetic_visualisation_{METHOD}.png", dpi=300)
plt.show()

# theoretical models for electric field and Coulomb potential.
logarithmic = lambda r, a, b: a * np.log(r) + b
inverse = lambda r, a, b: a / r + b

# fit curves with optimal parameters
opt_Az, _ = curve_fit(logarithmic, distance, Az_line)

# plot potential
plt.figure()
plt.plot(distance, Az_line, "o", label="Simulation", markersize=4)
plt.plot(distance, logarithmic(distance, *opt_Az), "-", label="Fit")
plt.xlabel(r"Distance from wire $r$")
plt.ylabel(r"Vector Potential $A_z$")
plt.title("Vector Potential vs. Distance")
plt.legend()
plt.savefig(f"./plots/vector_potential_{METHOD}.png", dpi=300)
plt.show()

# same as above, curve fitting
opt_B, _ = curve_fit(inverse, distance, B_field_line)

# plot E_field strength
plt.figure()
plt.plot(distance, B_field_line, "o", label="Simulation", markersize=4)
plt.plot(distance, inverse(distance, *opt_B), "-", label="Fit")
plt.xlabel(r"Distance from wire $r$")
plt.ylabel(r"Magnetic Field Strength $|\mathbf{B}|$")
plt.title("Magnetic Field Strength vs. Distance")
plt.legend()
plt.savefig(f"./plots/B_field_strength_{METHOD}.png", dpi=300)
plt.show()