# Import the champions
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

'''
LOAD DATA SET
'''

# Loading the data set
# Trim the trailing zero artifact
positions_raw = np.loadtxt('monomer_positions.txt', delimiter=',')
positions_array = positions_raw[:,0:-1]

# Determine how many different conformations there are
# Determine how many monomers there are
num_steps = positions_array.shape[0]
num_monomers = positions_array.shape[1] // 3

# Recast the positions array as a list of positions at every step
# Also, reshape the positions into a N x 3 array
shape = (num_monomers, 3)
positions_list = [np.reshape(row, shape) for row in positions_array]

'''
PLOT SYSTEM
'''

# Initiate axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Initial positions
graph = ax.scatter(positions_list[0][:,0], \
					positions_list[0][:,0], \
					positions_list[0][:,0],
					s=100, c='b')

# Update atom positions
def update(frame):

	graph._offsets3d = (positions_list[frame][:,0], \
						positions_list[frame][:,1], \
						positions_list[frame][:,2])

'''
CONFIGURE PLOT SETTINGS
'''

'''
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

'''
ANIMATE THE POLYMER
'''

fps = 60
ani = FuncAnimation(fig, update, num_steps, interval=1000/fps, blit=False)
plt.show()
