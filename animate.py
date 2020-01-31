# Import the champions
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

'''
LOAD DATA SET
'''

# Loading the data set
filename = 'data/monomer_positions.txt'
data = np.loadtxt(filename, delimiter=',')

# Isolate positions and end-to-end distances
positions_array = data[:,0:-1]
distances_array = data[:,-1]

# Determine number of conformations and monomers
num_steps = positions_array.shape[0]
num_monomers = positions_array.shape[1] // 3

# Recast the positions array as a list of positions at every step
# Also, reshape the positions into a N x 3 array
shape = (num_monomers, 3)
positions_list = [np.reshape(row, shape) for row in positions_array]

'''
PLOT SYSTEM
'''

# Initiate axes for animation and histogram
fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122)

# Initial positions
line = ax1.plot(positions_list[0][:,0], \
				positions_list[0][:,0], \
				positions_list[0][:,0],
				'-o')[0]
				
ax2.hist(distances_array[0], color='g', \
			bins=num_monomers, range=(0, num_monomers))
ax2.axvline(np.sqrt(num_monomers), zorder=1, \
			c='r', linewidth=3, label=r'R$_{RMS}$')

'''
https://matplotlib.org/gallery/animation/animated_histogram.html
'''

# Update axes
def update(frame):

	# Polymer conformation
	line.set_data(positions_list[frame][:,0:2].T)
	line.set_3d_properties(positions_list[frame][:,2])
	
	# Histogram
	ax2.hist(distances_array[0:frame+1], color='g', \
			bins=num_monomers, range=(0, num_monomers))

'''
CONFIGURE PLOT SETTINGS
'''

# Set limits for 3D plot
font = 16
limit = np.sqrt(num_monomers)
ax1.set_xlim3d([-limit, limit])
ax1.set_ylim3d([-limit, limit])
ax1.set_zlim3d([-limit, limit])
ax2.set_xlim([0, num_monomers])
ax2.legend(fontsize=font-4, loc='upper right')
ax1.set_title('Freely Jointed 3D Random Walk', fontsize=font)
ax2.set_title('End-to-end Distances', fontsize=font)

'''
ANIMATE THE POLYMER
'''

fps = 60
ani = FuncAnimation(fig, update, num_steps, interval=1000/fps, blit=False)
#ani.save('data/animation.gif', writer='imagemagick', fps=60)
plt.show()
