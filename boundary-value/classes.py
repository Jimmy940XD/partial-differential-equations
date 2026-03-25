import numpy as np
from functions import jacobi

class Lattice:
    def __init__(self, side_len: int):
        self.dx_ = 1
        self.dt_ = self.dx_ / 2
        self.lattice = np.zeros(shape=(side_len, side_len, side_len))

    def evolve(self):
        """
        Handler method. If assigned to a variable, it updates the
        lattice to the next iteration and returns absolute
        difference between the previous and the current lattice.
        If not assigned to a variable, it only updates the lattice.
        """
        return jacobi(self.lattice)