import numpy as np
from functions import jacobi

class Lattice:
    def __init__(self, side_len: int):
        self.dx_ = 1
        self.dt_ = self.dx_ / 2
        self.lattice = np.zeros(shape=(side_len, side_len, side_len))

    def evolve(self):
        jacobi(self.lattice)