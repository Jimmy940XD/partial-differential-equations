# partial-differential-equations
Checkpoint 3 of the Modelling and Visualisation in Physics course from The University of Edinburgh.

part 1.

optimal $\delta \tilde{t}$ values: at 0.03 and above, the algorithm diverges. Up to 0.025, the algorithm converges. The limiting point must be somewhere between 0.03 and 0.025.

$\phi _0$ values: negative and positive values show a predominance of oil and water (or the respective fluids used) concentration, respectively. They form droplets of the one that has least concentration. This two possible systems undergo ``Ostwald ripening``. A null initial value gives an equal concentration for both fluids, reaching a ratio of concentrations much more balanced, closer to 1. Worm-like shapes from both fluids can be seen, growing larger and thicker as they try to minimize the free energy between phases. This setup undergoes ``spinoidal decomposition``.

part 2.

specifying $\phi=0$ on the boundary corresponds to a Dirichlet BC, as we're specifying $\phi$ itself, not its derivative.

task 7: the two plots of potential and electric field vs. distance each decay as a power of an inverse the further they are from the centre of the lattice (where the point-charge is located). This agrees with what we would expect from Gauss's Law, which yields an inverse relation for the potential and an inverse-squared relation for the electric field.

task 8: the equation $\nabla^2\mathbf{A}=-\mu _0\mathbf{J}$ is actually a system of equations disguised as a single vector equation. The system reads:\
$$
\partial^2{A}{x}
$$