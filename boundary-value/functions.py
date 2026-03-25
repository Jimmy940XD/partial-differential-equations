import numpy as np
from numba import njit


@njit
def jacobi(lattice: np.ndarray):
    SIDE_LEN = lattice.shape[0]
    phi0 = lattice.copy()
    rho = np.zeros_like(lattice)
    phi = np.zeros_like(lattice)
    for i in range(SIDE_LEN):
        for j in range(SIDE_LEN):
            for k in range(SIDE_LEN):
                # calculate the next potential
                phi[i, j, k] = 1 / 6 * (phi[i + 1, j, k] + phi[i - 1, j, k] + \
                                        phi[i, j + 1, k] + phi[i, j - 1, k] + \
                                        phi[i, j, k + 1] + phi[i, j, k - 1] + \
                                        rho[i, j, k])
    lattice[:] = phi[:]
    return np.abs(phi - phi0)


@njit
def calc_Efield(phi: np.ndarray):
    SIDE_LEN = phi.shape[0]
    E_field = np.zeros_like(phi)
    for i in range(SIDE_LEN):
        for j in range(SIDE_LEN):
            for k in range(SIDE_LEN):
                # calculation of components of electric field using a centred difference
                Ex = phi[i + 1, j, k] - phi[i - 1, j, k]
                Ey = phi[i, j + 1, k] - phi[i, j - 1, k]
                Ez = phi[i, j, k + 1] - phi[i, j, k - 1]
                E_field[i, j, k] = -1 / 2 * np.array([Ex, Ey, Ez]) # the actual vector per cell
    return E_field