import numpy as np
from numba import njit, prange


@njit
def jacobi(info: tuple):
    """
    Updates the current lattice in-place using the Jacobi method.
    
    Returns
    -------
    max_diff: float
        The highest of the absolute differences between the previous
        step's lattice and the current step's lattice.
    """
    lattice, source = info
    SIDE_LEN = lattice.shape[0]
    lattice0 = lattice.copy()
    for i in range(1, SIDE_LEN - 1):
        for j in range(1, SIDE_LEN - 1):
            for k in range(1, SIDE_LEN - 1):
                # calculate the next potential
                lattice[i, j, k] = 1 / 6 * (lattice0[i + 1, j, k] + lattice0[i - 1, j, k] + \
                                        lattice0[i, j + 1, k] + lattice0[i, j - 1, k] + \
                                        lattice0[i, j, k + 1] + lattice0[i, j, k - 1] + \
                                        source[i, j, k])
    return np.max(np.abs(lattice - lattice0))


@njit
def gauss_seidel(info: tuple, over_relax: bool = False, omega: float = None):
    """
    Updates the current lattice in-place using the red-black Gauss-Seidel method.
    If ``over_relax`` is ``True``, applies SOR method.
    
    Returns
    -------
    max_diff: float
        The highest of the absolute differences between the previous
        step's lattice and the current step's lattice.
    """
    lattice, source = info
    SIDE_LEN = lattice.shape[0]
    lattice0 = lattice.copy()
    # the two loops allow for parallelisation using the red-black method
    for i in prange(1, SIDE_LEN - 1):
        for j in range(1, SIDE_LEN - 1):
            for k in range(1, SIDE_LEN - 1):
                red = (i + j + k) % 2 == 0
                if red:
                    lattice[i, j, k] = 1 / 6 * (lattice[i + 1, j, k] + lattice[i - 1, j, k] + \
                                            lattice[i, j + 1, k] + lattice[i, j - 1, k] + \
                                            lattice[i, j, k + 1] + lattice[i, j, k - 1] + \
                                            source[i, j, k])
                if over_relax: # apply SOR algorithm
                    lattice[i, j, k] = (1 - omega) * lattice0[i, j, k] + omega * lattice[i, j, k]

    for i in prange(1, SIDE_LEN - 1):
        for j in range(1, SIDE_LEN - 1):
            for k in range(1, SIDE_LEN - 1):
                red = (i + j + k) % 2 == 0
                if not red:
                    lattice[i, j, k] = 1 / 6 * (lattice[i + 1, j, k] + lattice[i - 1, j, k] + \
                                            lattice[i, j + 1, k] + lattice[i, j - 1, k] + \
                                            lattice[i, j, k + 1] + lattice[i, j, k - 1] + \
                                            source[i, j, k])
                    if over_relax: # apply SOR algorithm
                        lattice[i, j, k] = (1 - omega) * lattice0[i, j, k] + omega * lattice[i, j, k]
    return np.max(np.abs(lattice - lattice0))


@njit
def sor(info: tuple):
    """
    Updates the current lattice in-place using the SOR method.
    """
    lattice, source, omega = info
    lattice0 = lattice.copy()
    info = (lattice, source) # redefine info for better readability
    gauss_seidel(info, over_relax=True, omega=omega)
    return np.max(np.abs(lattice - lattice0))


@njit
def calc_Efield(phi: np.ndarray):
    """
    Calculates and returns the 3D electric field at each cell of
    the lattice based on the current step's Coulomb potential.
    """
    SIDE_LEN = phi.shape[0]
    E_field = np.zeros(shape=(SIDE_LEN, SIDE_LEN, SIDE_LEN, 3)) # the 3 is for the three-dimensional vector in each cell
    for i in range(1, SIDE_LEN - 1):
        for j in range(1, SIDE_LEN - 1):
            for k in range(1, SIDE_LEN - 1):
                # calculation of components of electric field using a centred difference
                # the i and j indices are swapped because i represents y and j represents x in code
                Ex = phi[i, j + 1, k] - phi[i, j - 1, k]
                Ey = phi[i + 1, j, k] - phi[i - 1, j, k]
                Ez = phi[i, j, k + 1] - phi[i, j, k - 1]
                # the actual vector per cell
                E_field[i, j, k, 0] = -1 / 2 * Ex
                E_field[i, j, k, 1] = -1 / 2 * Ey
                E_field[i, j, k, 2] = -1 / 2 * Ez
    return E_field


@njit
def calc_Bfield(Az: np.ndarray):
    """
    Calculates and returns the 3D magnetic field at each cell of
    the lattice based on the current step's vector potential.
    """
    SIDE_LEN = Az.shape[0]
    B_field = np.zeros(shape=(SIDE_LEN, SIDE_LEN, SIDE_LEN, 2))
    for i in range(1, SIDE_LEN - 1):
        for j in range(1, SIDE_LEN - 1):
            for k in range(1, SIDE_LEN - 1):
                # calculation of components of magnetic field using a centred difference
                # the i and j indices aren't swapped because i represents y and j represents x in code
                Bx = Az[i + 1, j, k] - Az[i - 1, j, k]
                By = -(Az[i, j + 1, k] - Az[i, j - 1, k])
                # the actual vector per cell
                B_field[i, j, k, 0] = 1 / 2 * Bx
                B_field[i, j, k, 1] = 1 / 2 * By
    return B_field