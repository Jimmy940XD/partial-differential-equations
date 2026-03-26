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

distances = []
potentials = []
E_field_strengths = []
convergence = False
while not convergence:
    distance = state.evolve(METHOD) # read method's docstring for deeper understanding
    distances.append(distance)
    potential = np.max(np.abs(state.lattice))
    potentials.append(potential)
    E_field_squared = state.E_field**2 # square E_field components
    E_field_strength = np.sqrt(np.sum(E_field_squared, axis=3)) # takes magnitude of E_field
    E_field_strengths.append(E_field_strength)
    if distance <= TOLERANCE:
        convergence = True

MID = SIDE_LEN // 2
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
plt.plot(distances, potentials)
plt.xlabel("Distance")
plt.ylabel(r"Potential $\phi$")
plt.title("Potential vs. Distance")
plt.savefig(f"./plots/potential_{METHOD}.png", dpi=300)
plt.show()

# plot E_field strength
plt.figure()
plt.plot(distances, E_field_strengths)
plt.xlabel("Distance")
plt.ylabel(r"Electric Field Strength $|\underline{E}|$")
plt.title("Electric Field Strength vs. Distance")
plt.savefig(f"./plots/E_field_strength_{METHOD}.png", dpi=300)
plt.show()