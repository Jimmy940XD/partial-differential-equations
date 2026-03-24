import numpy as np
from functions import cahn_hilliard_cfd

class Lattice:
    def __init__(self, side_len: int, phi0: float):
        # set-up of the initial lattice
        base = np.zeros(shape=(side_len, side_len)).fill(phi0)
        noise = np.random.uniform(-.01, .01, size=(side_len, side_len))
        self.lattice = base + noise

    def evolve(self):
        cahn_hilliard_cfd(self.lattice)