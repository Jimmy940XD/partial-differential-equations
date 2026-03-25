from classes import Lattice
import matplotlib.pyplot as plt

SIDE_LEN = 100
PHI0 = float(input("Type an initial value for phi_0: "))

state = Lattice(SIDE_LEN, PHI0)

energies = []
times = []
n = 0
convergence = False
while not convergence:
    times.append(n)
    energies.append(state.free_energy)
    state.evolve()
    if n >= 50000:
        convergence = True
    n += 1

plt.plot(times, energies)
plt.xlabel(r"Timestep $t^{\tilde}$")
plt.ylabel(r"Free energy $f^{\tilde}")
plt.title("Free energy versus time")
plt.show()