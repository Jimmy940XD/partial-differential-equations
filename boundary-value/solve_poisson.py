from classes import Lattice
import numpy as np
import matplotlib.pyplot as plt
import os

SIDE_LEN = int(input("Choose the system's side length: "))
METHOD = input("Choose the evolution method (Jacobi/Gauss-Seidel): ").lower()
TOLERANCE = float(input("Choose the accuracy of the solution: "))
# choose rho as one likes
rho = None # this tells the state to initialize with a monopole configuration

state = Lattice(SIDE_LEN, rho)

convergence = False
while not convergence:
    diff = state.evolve(METHOD) # read method's docstring for deeper understanding
    if diff <= TOLERANCE:
        convergence = True

MID = SIDE_LEN // 2
E_field_squared = state.E_field**2 # square E_field components
E_field_strength = np.sqrt(np.sum(E_field_squared, axis=3)) # takes magnitude of E_field
distance = np.arange(1, -2) # distance from cell next to CENTER up to boundary (excluding it)
# lines depndent only on distance to charge
phi_line = state.lattice[MID + 1:-1, MID, MID]
E_field_line = E_field_strength[MID + 1:-1, MID, MID]

# get data at midplane cut
phi_midplane = state.lattice[:, :, MID]
Ex_midplane = state.E_field[:, :, MID, 0]
Ey_midplane = state.E_field[:, :, MID, 1]

# save data
os.makedirs("./datafiles", exist_ok=True)
np.savetxt("./datafiles/potential.dat", phi_midplane, header="Potential Lattice")
np.savetxt("./datafiles/E_field.dat", state.E_field[:, :, MID].reshape(-1, 3), header="Electric Field Lattice")

x = np.arange(SIDE_LEN)
y = np.arange(SIDE_LEN)
X, Y = np.meshgrid(x, y)
fig, ax = plt.subplots()
contour = ax.contourf(X, Y, phi_midplane, levels=50, cmap="RdBu_r") # creates the contour with 50 discretized levels
colbar = plt.colorbar(contour, ax=ax, fraction=.046, pad=.04)
colbar.set_label(r"Coulomb potential ($\varphi$)")
jump = 4 # reduces the number of arrows that will be seen in the E_field plot
ax.quiver(X[::jump, ::jump], Y[::jump, ::jump],
        Ex_midplane[::jump, ::jump], Ey_midplane[::jump, ::jump],
        color="white",
        pivot="mid",
        alpha=.8) # creates a 2D slice of the 3D E_field
ax.set_title("Cut through midplane: Potential and Electric Field")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal") # ensures the image doesn't get stretched
plt.show()

# plot potential
plt.figure()
plt.plot(distance, phi_line)
plt.xlabel(r"Distance from charge $r$")
plt.ylabel(r"Potential $\varphi$")
plt.title("Potential vs. Distance")
plt.savefig(f"./plots/potential_{METHOD}.png", dpi=300)
plt.show()

# plot E_field strength
plt.figure()
plt.plot(distance, E_field_line)
plt.xlabel(r"Distance from charge $r$")
plt.ylabel(r"Electric Field Strength $|\mathbf{E}|$")
plt.title("Electric Field Strength vs. Distance")
plt.savefig(f"./plots/E_field_strength_{METHOD}.png", dpi=300)
plt.show()