import numpy as np
import functions as f

class Lattice:
    """
    General 3D lattice.
    
    Attributes
    ----------
    self.lattice: numpy.ndarray
        3D array representing a lattice.
    self.E_field: numpy.ndarray
        Property method. Returns an array of 3D electric field vectors
        at each cell of the lattice.
    """
    def __init__(self, side_len: int, source: np.ndarray = None):
        self.lattice = np.zeros(shape=(side_len, side_len, side_len))
        if source is None:
            # make it point charge by default
            self.source = np.zeros_like(self.lattice)
            MID = side_len // 2
            CENTER = (MID, MID, MID)
            self.source[CENTER] = 1
        else:
            self.source = source

    def evolve(self, method: str):
        """
        Handler method. If assigned to a variable, it updates the
        lattice to the next iteration and returns the highest absolute
        difference among all the elements between the previous and the
        current lattice. If not assigned to a variable, it only updates
        the lattice.
        """
        info = (self.lattice, self.source)
        method = method.lower()
        if method == "jacobi":
            return f.jacobi(info)
        if method == "gauss-seidel":
            return f.gauss_seidel(info)
    
    @property
    def E_field(self):
        """
        Handler method. Returns an array of electric field 3D vectors at
        each cell of the lattice.
        """
        return f.calc_Efield(self.lattice)