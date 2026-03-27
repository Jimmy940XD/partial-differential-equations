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
    self.B_field: numpy.ndarray
        Property method. Returns an array of 3D magnetic field vectors
        at each cell of the lattice.
    """
    def __init__(self, side_len: int, source: np.ndarray):
        self.source = source
        self.lattice = np.zeros(shape=(side_len, side_len, side_len))

    def evolve(self, omega: float, method: str = "sor"):
        """
        Handler method. If assigned to a variable, it updates the
        lattice to the next iteration and returns the highest absolute
        difference among all the elements between the previous and the
        current lattice. If not assigned to a variable, it only updates
        the lattice.

        Parameters
        ----------
        method: str, optional
            The evolution algorithm to be used. By default, it is set to
            the SOR method, which is the fastest.
        omega: float, optional
            The omega parameter for the evolution of the system using SOR.
            Only needed for SOR, hence it is ``None`` by default as it is
            expected that if no value of ``omega`` is passed, ``method``
            will be. 
        """
        if omega is None:
            info = (self.lattice, self.source)
            method = method.lower()
            if method == "jacobi":
                return f.jacobi(info)
            if method == "gauss-seidel":
                return f.gauss_seidel(info)
        else:
            info = (self.lattice, self.source, omega)
            return f.sor(info)
    
    @property
    def E_field(self):
        """
        Handler method. Returns an array of electric field 3D vectors at
        each cell of the lattice.
        """
        return f.calc_Efield(self.lattice)
    
    @property
    def B_field(self):
        """
        Handler method. Returns an array of magnetic field 3D vectors at
        each cell of the lattice.
        """
        return f.calc_Bfield(self.lattice)