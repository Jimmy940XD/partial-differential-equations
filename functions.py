import numpy as np
from numba import njit


@njit
def cahn_hilliard_cfd(info: tuple):
    """
    Updating function following the Cahn-Hilliard equation, using
    a centred spatial difference and a forward time difference.
    """
    lattice, DX_, dt_ = info
    SIDE_LEN = lattice.shape[0]
    phi = lattice.copy()
    mu = np.zeros_like(lattice)
    for i in range(SIDE_LEN):
        for j in range(SIDE_LEN):
            mu[i, j] = -phi[i, j] * (1 - phi[i, j]**2) - \
                        1 / DX_**2 * (phi[(i + 1) % SIDE_LEN, j] + \
                        phi[(i - 1) % SIDE_LEN, j] + phi[i, (j + 1) % SIDE_LEN] + \
                        phi[i, (j - 1) % SIDE_LEN] - 4 * phi[i, j])
    for i in range(SIDE_LEN):
        for j in range(SIDE_LEN):
            phi[i, j] = phi[i, j] + dt_ / DX_**2 * (mu[(i + 1) % SIDE_LEN, j] + \
                        mu[(i - 1) % SIDE_LEN, j] + mu[i, (j + 1) % SIDE_LEN] + \
                        mu[i, (j - 1) % SIDE_LEN] - 4 * mu[i, j])
    lattice[:] = phi[:]


@njit
def free_energy(info: tuple):
    """
    Calculates the dimesnionless free energy density.
    """
    phi, DX_ = info
    SIDE_LEN = phi.shape[0]
    f = np.zeros_like(phi)
    for i in range(SIDE_LEN):
        for j in range(SIDE_LEN):
            f[i, j] = -1 / 2 * phi[i, j]**2 + 1 / 4 * phi[i, j]**4 + \
                        1 / 2 * (((phi[(i + 1) % SIDE_LEN, j] - \
                        phi[(i - 1) % SIDE_LEN, j]) / (2 * DX_))**2 + \
                        ((phi[i, (j + 1) % SIDE_LEN] - phi[i, (j - 1) % SIDE_LEN]) / (2 * DX_))**2)
    return f