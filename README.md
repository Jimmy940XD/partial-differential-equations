# Partial Differential Equations in Boundary and Initial Value problems
**Checkpoint 3 - Modelling and Visualisation in Physics**

## Description
Checkpoint 3 of the Modelling and Visualisation in Physics course from The University of Edinburgh.
This project explores several algorithms that discretize partial differential equations in order to solve problems with boundary conditions or initial conditions.

## Libraries needed
* `NumPy`
* `SciPy`
* `Numba`
* `os`
* `Matplotlib`

## Files and Directories
* `initial-value`: Folder for the initial value problems. Contains:
    * `classes.py`: This file is intended for definition of the `Lattice` class that will be used in all initial value problems.
    * `functions.py`: This file stores all the functions and algorithms that will be used in these problems. It uses a Just-In-Time compiler to speed up the algorithms.
    * `animation.py`: This file runs a program that shows an animation of the evolution of a fluid-mixing system following an initial condition.
    * `free_energy.py`: This file plots the free energy of the system at each timestep to show how it converges when $t\rightarrow\infty$.
* `boundary-value`: Folder for the boundary value problems. Contains:
    * `classes.py`: Same as in `initial-value`. Even though the class is called the same, it is a complete different object in memory.
    * `functions.py`: Same as in `initial-value`.
    * `solve_poisson.py`: This file solves the Poisson equation for a user input of the system's side length, the algorithm to be used to solve it, and the accuracy of the solution (the tolerance to convergence). It plots a contour plot of the Coulomb potential for a specific charge distribution (to be modified as desired at the beginning of the script) and a vector plot of the electric field. It also plots two sampled and fitted curves that show the strength of the potential and field against distance to the charge distribution. A datafile with the values for the field and potential is created as well.
    * `solve_magnetic.py`: Same as previous file, but with the magnetic field and magnetic vector potential.
    * `convergence.py`: Contains a program that shows convergence of the three algorithms available, and how fast each of them converges (how many iterations are needed).
    * `minimisation.py`: File containing a program that determines the value for $\omega$ that minimizes the time it takes the SOR algorithm to converge. It also shows a plot of number of iterations agains $\omega$ values.

## How to Run the Code
- Make sure you have a `plots/` directory created in each of the `...-value` folders. Otherwise the plots will not be saved and a ``FileNotFoundError`` will be raised by Python. This is the only data issue to be concerned about.
- To get the plots and datafiles, just run the corresponding scripts.
- In `minimisation.py`, if a larger than 50 side length is chosen the script will run for several minutes up to hours (depending on how big it is). The recommended size for optimal results in short time is between 25 and 50. It should only take up to 5 minutes at worst. This will also happen if the ``STEP`` variable is hard-coded with a very small value (the current hard-coded is optimal).
- For the animation, just run the `animation.py` script.

## Tasks

### Initial Value Problems
* Task 3: Setting $\delta\tilde{x} = 1$ for cell discretization I found some optimal $\delta\tilde{t}$ values: at 0.03 and above, the algorithm diverges. Up to 0.025, the algorithm converges. The limiting point must be somewhere between 0.03 and 0.025.
* Task 4: Exploring $\varphi _0$ values one finds some qualitative behaviours. Negative and positive values show a predominance of oil and water (or the respective fluids used) concentration, respectively. They form droplets of the one that has least concentration. This two possible systems undergo ``Ostwald ripening``. A null initial value gives an equal concentration for both fluids, reaching a ratio of concentrations much more balanced, closer to 1. Worm-like shapes from both fluids can be seen, growing larger and thicker as they try to minimize the free energy between phases. This setup undergoes ``spinoidal decomposition``.

### Boundary Value Problems
* Task 1: specifying $\varphi=0$ on the boundary corresponds to a Dirichlet BC, as we're specifying $\varphi$ itself, not its derivative.
* Task 7: the two plots of potential and electric field vs. distance each decay as an inverse law the further they are from the centre of the lattice (where the point-charge is located). This agrees with what we would expect from Gauss's Law, which yields an inverse relation for the potential and an inverse-squared relation for the electric field.
* Task 8: the equation $\nabla^2\mathbf{A}=-\mu _0\mathbf{J}$ is actually a system of equations disguised as a single vector equation. The system reads:

$$
\begin{align*}
\nabla^2A_x &= -\mu _0J_x, \\
\nabla^2A_y &= -\mu _0J_y, \\
\nabla^2A_z &= -\mu _0J_z.
\end{align*}
$$

Since $\mathbf{J}=(0, 0, J_z)$, the RHS of the first two equations becomes zero instantly, implying $A_{x,y}=0$ because the solution of the Laplace equation is that in a physical space without boundaries (or boundaries set to zero, like in our case).

Hence we have

$$
\frac{\partial^2 A_z}{\partial x^2}+\frac{\partial^2 A_z}{\partial y^2}+
\frac{\partial^2 A_z}{\partial z^2}=-\mu _0 J_z.
$$

Since the current only goes through the $z$-direction, the system looks the same no matter what $z$ coordinate you are in. Therefore the derivative of $A_z$ with respect to $z$ is zero, as $J_z$ is independent of $z$. As a result, $A_z$ is only $(x,y)$-dependent.

Good boundary conditions, as was mentioned above, are to take the boundaries equal to zero, as was done with Poisson's equation and the Coulomb potential. These are Dirichlet BCs.