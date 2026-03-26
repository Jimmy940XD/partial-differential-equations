from classes import Lattice

SIDE_LEN = 100
TOLERANCE = 1e-5

state = Lattice(SIDE_LEN)

METHODS = ("Jacobi", "Gauss-Seidel")
for method in METHODS:
    convergence = False
    while not convergence:
        distance = state.evolve(method) # read method's docstring for deeper understanding
        if distance <= TOLERANCE:
            convergence = True
    print(f"Convergence reached using {method} at distance: {distance}")