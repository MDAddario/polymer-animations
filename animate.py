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
filename = 'data/monomer_positions.txt'
positions_raw = np.loadtxt(filename, delimiter=',')
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
line = ax.plot(positions_list[0][:,0], \
				positions_list[0][:,0], \
				positions_list[0][:,0],
				'-o')[0]

# Update atom positions
def update(frame):

	line.set_data(positions_list[frame][:,0:2].T)
	line.set_3d_properties(positions_list[frame][:,2])

'''
CONFIGURE PLOT SETTINGS
'''

# Set axis limits
limit = np.sqrt(num_monomers) / 2
ax.set_xlim3d([-limit, limit])
ax.set_ylim3d([-limit, limit])
ax.set_zlim3d([-limit, limit])
ax.set_title('Freely Jointed 3D Random Walk', fontsize=20)

'''
ANIMATE THE POLYMER
'''

fps = 60
ani = FuncAnimation(fig, update, num_steps, interval=1000/fps, blit=False)
plt.show()
