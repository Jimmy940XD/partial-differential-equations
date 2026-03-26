from classes import Lattice

SIDE_LEN = 100
TOLERANCE = 1e-5

# this loop initializes and runs the system for each method
METHODS = ("Jacobi", "Gauss-Seidel")
for method in METHODS:
    state = Lattice(SIDE_LEN)
    n = 0
    convergence = False
    while not convergence:
        distance = state.evolve(method) # read method's docstring for deeper understanding
        n += 1
        if distance <= TOLERANCE:
            convergence = True
    print(f"{method} converged after {n} iterations at distance: {distance}")