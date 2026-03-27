from classes import Lattice
import numpy as np
import matplotlib.pyplot as plt

SIDE_LEN = int(input("Choose the side length of the system (recommended between 25 and 50): "))
TOLERANCE = float(input("Choose the accuracy of the solution: "))
STEP = .01 # this determines how many omega values are checked. Keep in mind the smaller the more calculation time
# current is chosen because it takes longer than point charge to converge with Gauss-Seidel
# main goal is to see clear faster convergence but edit as desired
current = np.zeros(shape=(SIDE_LEN, SIDE_LEN, SIDE_LEN))
MID = SIDE_LEN // 2
current[MID, MID, 1:-1] = 1

iterations = []
omegas = np.concatenate((np.arange(1, 1.7, .1), np.arange(1.7, 1.99, STEP))) # initialize omega values
for omega in omegas:
    state = Lattice(SIDE_LEN, current)
    convergence = False
    n = 0
    while not convergence:
        diff = state.evolve(omega=omega)
        n += 1
        if diff <= TOLERANCE:
            convergence = True
    iterations.append(n)
index_min = np.argmin(iterations) # get the index for the minimum value
omega_min = omegas[index_min]
print(f"With 1 <= ω <= 2 and {STEP} separation, ω = {round(omega_min, 2)} shortens convergence time the most.")

plt.figure()
plt.plot(omegas, iterations)
plt.xlabel(r"$\omega$")
plt.ylabel(r"# iteration $n$")
plt.title(r"# iterations per each value of $\omega$")
plt.savefig(f"./plots/omega.png", dpi=300)
plt.show()