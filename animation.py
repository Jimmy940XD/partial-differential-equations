import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import ListedColormap
from classes import Lattice

def update(frame):
    state.evolve()

    img.set_data(state.lattice)
    return [img]

SIDE_LEN = 50
PHI0 = float(input("Type an initial value for phi_0: "))

state = Lattice(SIDE_LEN, PHI0)

fig, ax = plt.subplots()
cmap = ListedColormap(["yellow", "blue"])
img = ax.imshow(state.lattice, cmap=cmap, interpolation="none", vmin=(PHI0 - .01), vmax=(PHI0 + .01))

animation = FuncAnimation(fig, update, frames=None, interval=50, blit=True, cache_frame_data=False)
plt.show()