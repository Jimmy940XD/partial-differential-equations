from classes import Lattice
import numpy as np

SIDE_LEN = int(input("Choose the side length of the system: "))
TOLERANCE = float(input("Choose the accuracy of the solution: "))
STEP = .01
# current is chosen because it takes longer than point charge to converge with Gauss-Seidel
# main goal is to see clear faster convergence but edit as desired
current = np.zeros(shape=(SIDE_LEN, SIDE_LEN, SIDE_LEN))
MID = SIDE_LEN // 2
current[MID, MID, 1:-1] = 1

state = Lattice(SIDE_LEN, current)

iterations = []
omegas = np.arange(1, 2 + STEP, STEP) # initialize omega values
for omega in omegas:
    convergence = False
    n = 0
    while not convergence:
        diff = state.evolve(omega)
        n += 1
        if diff <= TOLERANCE:
            convergence = True
    iterations.append(n)
index_min = np.argmin(iterations) # get the index for the minimum value
omega_min = omegas[index_min]
print(f"With 1 <= ω <= 2 and {STEP} separation, ω = {omega_min} shortens convergence time the most.")