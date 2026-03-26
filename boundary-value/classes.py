import numpy as np
import functions as f

class Lattice:
    def __init__(self, side_len: int):
        self.lattice = np.zeros(shape=(side_len, side_len, side_len))

    def evolve(self, method: str):
        """
        Handler method. If assigned to a variable, it updates the
        lattice to the next iteration and returns the highest absolute
        difference among all the elements between the previous and the
        current lattice. If not assigned to a variable, it only updates
        the lattice.
        """
        method = method.lower()
        if method == "jacobi":
            return f.jacobi(self.lattice)
        if method == "gauss-seidel":
            return f.gauss_seidel(self.lattice)
    
    @property
    def E_field(self):
        """
        Handler method. Returns an array of electric field 3D vectors at
        each cell of the lattice.
        """
        return f.calc_Efield(self.lattice)