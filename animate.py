# Import the champions
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

'''
Plotting the system
'''

# Initiate axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Initial positions
graph = ax.scatter(positions[num_atoms * 0: num_atoms * 1, 0], \
					positions[num_atoms * 1: num_atoms * 2, 0], \
					positions[num_atoms * 2: num_atoms * 3, 0], \
					s=100, c='b')

# Update atom positions
def update(frame):

	graph._offsets3d = (positions[num_atoms * 0: num_atoms * 1, frame], \
						positions[num_atoms * 1: num_atoms * 2, frame], \
						positions[num_atoms * 2: num_atoms * 3, frame])

# Set axis limits
ax.set_xlim3d([-a, a * N_1])
ax.set_ylim3d([-b, b * N_2])
ax.set_zlim3d([-c, c * N_3])
ax.set_title('Cubic Lattice Phonons', fontsize=20, y=1.08)

# Hide grid lines and ticks
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

'''
Animate!
'''

ani = FuncAnimation(fig, update, tot_frames, interval=1000/fps, blit=False)
plt.show()
