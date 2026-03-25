import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap
from classes import Lattice

def update(frame):
    for _ in range(100):
        state.evolve()

    img.set_data(state.lattice)
    return [img]

SIDE_LEN = 100
PHI0 = float(input("Type an initial value for phi_0: "))

state = Lattice(SIDE_LEN, PHI0)

fig, ax = plt.subplots()
cmap = LinearSegmentedColormap.from_list("yellow_blue", ["yellow", "blue"])
img = ax.imshow(state.lattice, cmap=cmap, interpolation="bilinear", vmin=-1, vmax=1)

animation = FuncAnimation(fig, update, frames=None, interval=1, blit=True, cache_frame_data=False)
plt.show()