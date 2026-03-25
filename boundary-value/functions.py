import numpy as np
from numba import njit


@njit
def jacobi(phi: np.ndarray):
    SIDE_LEN = phi.shape[0]
    phi0 = phi.copy()
    rho = np.zeros_like(phi)
    mid = SIDE_LEN // 2
    centre = (mid, mid, mid)
    rho[centre] = 1 # initialize as one likes
    for i in range(1, SIDE_LEN - 1):
        for j in range(1, SIDE_LEN - 1):
            for k in range(1, SIDE_LEN - 1):
                # calculate the next potential
                phi[i, j, k] = 1 / 6 * (phi[i + 1, j, k] + phi[i - 1, j, k] + \
                                        phi[i, j + 1, k] + phi[i, j - 1, k] + \
                                        phi[i, j, k + 1] + phi[i, j, k - 1] + \
                                        rho[i, j, k])
    return np.max(np.abs(phi - phi0))


@njit
def calc_Efield(phi: np.ndarray):
    SIDE_LEN = phi.shape[0]
    E_field = np.zeros((SIDE_LEN, SIDE_LEN, SIDE_LEN, 3)) # the 3 is for the three-dimensional vector in each cell
    for i in range(1, SIDE_LEN - 1):
        for j in range(1, SIDE_LEN - 1):
            for k in range(1, SIDE_LEN - 1):
                # calculation of components of electric field using a centred difference
                Ex = phi[i + 1, j, k] - phi[i - 1, j, k]
                Ey = phi[i, j + 1, k] - phi[i, j - 1, k]
                Ez = phi[i, j, k + 1] - phi[i, j, k - 1]
                # the actual vector per cell
                E_field[i, j, k, 0] = -1 / 2 * Ex
                E_field[i, j, k, 1] = -1 / 2 * Ey
                E_field[i, j, k, 2] = -1 / 2 * Ez
    return E_field