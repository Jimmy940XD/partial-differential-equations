import numpy as np
from functions import cahn_hilliard_cfd, free_energy

class Lattice:
    """
    A representation of a mixed-fluid lattice.

    Properties
    ----------
    self.dx_: int
        The discretized spatial step. It is hard-coded to 1 for
        simplicity (a cell is the discretisation).
    self.dt_: float
        The discretized timestep.
    self.lattice: numpy.ndarray
        The lattice.
    self.free_energy: numpy.ndarray
        Property method. Returns an array with the value of the free energy
        at each point in space.
    """
    def __init__(self, side_len: int, phi0: float):
        # discretization
        self.dx_ = 1
        self.dt_ = float(input("Enter the time discretisation: "))
        # set-up of the initial lattice
        base = np.zeros(shape=(side_len, side_len))
        base.fill(phi0)
        noise = np.random.uniform(-.01, .01, size=(side_len, side_len))
        self.lattice = base + noise


    def evolve(self):
        info = (self.lattice, self.dx_, self.dt_)
        cahn_hilliard_cfd(info)

    
    @property
    def free_energy(self):
        info = self.lattice, self.dx_
        return free_energy(info)