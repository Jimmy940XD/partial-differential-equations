from classes import Lattice

SIDE_LEN = 100

state = Lattice(SIDE_LEN)

convergence = False
while not convergence:
    distance = state.evolve("jacobi") # read method's docstring for deeper understanding
    tolerance = 1e-5
    if distance <= tolerance:
        convergence = True
print(f"Convergence reached using Jacobi at distance: {distance}")