# Import the champions
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider, Button, RadioButtons
import time
from numba import njit

'''
User set parameters
'''

# Number of unit cells along each direction
N_1 = 6
N_2 = 6
N_3 = 6

# Animation parameters
fps   = 60.0
T_max = 1.0

'''
Determine total number of frames in the animation (cubic lattice)
'''

# Range of possible reciprocal coefficients
m_1_range = np.arange(N_1)
m_2_range = np.arange(N_2)
m_3_range = np.arange(N_3)

# List of phonon frequencies
omega_1 = 2 * np.sin(m_1_range / N_1 * np.pi)
omega_2 = 2 * np.sin(m_2_range / N_2 * np.pi)
omega_3 = 2 * np.sin(m_3_range / N_3 * np.pi)

# Maximum phonon frequencies
omega_1_max = np.max(omega_1)
omega_2_max = np.max(omega_2)
omega_3_max = np.max(omega_3)

# Number of frames required for a full oscillation of given mode
tot_frames_1 = np.rint(fps / omega_1[1:] * omega_1_max * T_max).astype('int')
tot_frames_2 = np.rint(fps / omega_2[1:] * omega_2_max * T_max).astype('int')
tot_frames_3 = np.rint(fps / omega_3[1:] * omega_3_max * T_max).astype('int')

# Lowest common multiple for given branch
lcm_frames_1 = np.lcm.reduce(tot_frames_1)
lcm_frames_2 = np.lcm.reduce(tot_frames_2)
lcm_frames_3 = np.lcm.reduce(tot_frames_3)

# Lowest common multiple across all branches and phonons
tot_frames = np.lcm.reduce([lcm_frames_1, lcm_frames_2, lcm_frames_3])

'''
Position array layout:

frame = 0	frame = 1	...
[
x_atom_1	x_atom_1	...
...
x_atom_n	x_atom_n	...
y_atom_1	y_atom_1	...
...
y_atom_n	y_atom_n	...
z_atom_1	z_atom_1	...
...
z_atom_n	z_atom_n	...
]

'''

# Declare position array
num_atoms = N_1 * N_2 * N_3
positions = np.empty((3 * num_atoms, tot_frames))

# Enter atomic coordinates into position array
@njit
def set_atomic_positions(positions, index, x_data, y_data, z_data):
	
	positions[index + 0 * num_atoms, :] = x_data
	positions[index + 1 * num_atoms, :] = y_data
	positions[index + 2 * num_atoms, :] = z_data
	
# Add atomic motions to position array
@njit
def add_atomic_positions(positions, index, x_data, y_data, z_data):
	
	positions[index + 0 * num_atoms, :] += x_data
	positions[index + 1 * num_atoms, :] += y_data
	positions[index + 2 * num_atoms, :] += z_data

'''
Construct crystal lattice (currently fixed as cubic system)
'''

# Lattice spacing along each direction
a = 1.0
b = 1.0
c = 1.0

'''
Precompute indices
'''
indices = np.empty((N_1, N_2, N_3), dtype='int')

for n_1 in range(N_1):
	for n_2 in range(N_2):
		for n_3 in range(N_3):

			indices[n_1, n_2, n_3] = n_1 + n_2 * N_1 + n_3 * (N_1 * N_2)

# Place each atom in its equilibrium position
@njit
def set_lattice_offsets(positions):

	for n_1 in range(N_1):
		for n_2 in range(N_2):
			for n_3 in range(N_3):

				index = indices[n_1, n_2, n_3]
				set_atomic_positions(positions, index, a * n_1, b * n_2, c * n_3)
				
'''
PRECOMPUTE THE SPATIAL ARRAY COMPONENTS
'''
spatial_n1 = np.empty((num_atoms, 1))
spatial_n2 = np.empty((num_atoms, 1))
spatial_n3 = np.empty((num_atoms, 1))

for n_1 in range(N_1):
	for n_2 in range(N_2):
		for n_3 in range(N_3):
			
			index = indices[n_1, n_2, n_3]
			spatial_n1[index, 0] = n_1 / N_1
			spatial_n2[index, 0] = n_2 / N_2
			spatial_n3[index, 0] = n_3 / N_3

# Excite a photon within the system
@njit
def excite_phonon(positions, branch, covar_coeffs, amplitude):

	# Extract reciprocal coefficients
	m_1, m_2, m_3 = covar_coeffs
	
	# Select branch evec and eval
	if branch == 1:
		epsilon   = amplitude * np.array([1, 0, 0])
		omega     = omega_1[m_1]
		omega_max = omega_1_max
	
	elif branch == 2:
		epsilon   = amplitude * np.array([0, 1, 0])
		omega     = omega_2[m_2]
		omega_max = omega_2_max
	
	elif branch == 3:
		epsilon   = amplitude * np.array([0, 0, 1])
		omega     = omega_3[m_3]
		omega_max = omega_3_max
	
	# Configure temporal portion of vibration
	temporal = omega / fps / omega_max / T_max * np.arange(tot_frames)
	
	# Configure spatial portion of vibration
	spatial = m_1 * spatial_n1 + m_2 * spatial_n2 + m_3 * spatial_n3
	
	# Compute vibration and update positions
	oscillation = np.sin(2 * np.pi * (spatial - temporal))
	x_disp = epsilon[0] * oscillation
	y_disp = epsilon[1] * oscillation
	z_disp = epsilon[1] * oscillation
	
	positions += np.vstack((x_disp, y_disp, z_disp))

set_lattice_offsets(positions)

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
Slider fun!
'''

# Radio variable for branch selection
radio_branch = 1

# Setup the plot
plt.subplots_adjust(left=0.25, bottom=0.27)
ax.margins(x=0)
axcolor_1 = 'turquoise'
axcolor_2 = 'lightgoldenrodyellow'
axcolor_3 = 'lavender'

# Reciprocal lattice vector covariant coefficient sliders (left, bottom, width, height)
ax_amp = plt.axes([0.25, 0.22, 0.65, 0.03], facecolor=axcolor_2)
ax_m_1 = plt.axes([0.25, 0.17, 0.65, 0.03], facecolor=axcolor_2)
ax_m_2 = plt.axes([0.25, 0.12, 0.65, 0.03], facecolor=axcolor_2)
ax_m_3 = plt.axes([0.25, 0.07, 0.65, 0.03], facecolor=axcolor_2)

slid_amp = Slider(ax_amp, 'Amplitude', 0.0, 0.5, valinit=0.3, valstep=0.01)
slid_m_1 = Slider(ax_m_1, r'$m_1$', 0, N_1-1, valinit=0, valstep=1)
slid_m_2 = Slider(ax_m_2, r'$m_2$', 0, N_2-1, valinit=0, valstep=1)
slid_m_3 = Slider(ax_m_3, r'$m_3$', 0, N_3-1, valinit=0, valstep=1)

def update_phonon(val):
	
	# Retrieve slider values
	amp = slid_amp.val
	m_1 = int(slid_m_1.val)
	m_2 = int(slid_m_2.val)
	m_3 = int(slid_m_3.val)

	set_lattice_offsets(positions)

	# Benchmark phonon
	global radio_branch
	start = time.time()
	excite_phonon(positions, radio_branch, [m_1, m_2, m_3], amp)
	end = time.time()
	print('Excitation took: {:.2e}s'.format(end - start))
	
	# Update title
	if m_1 == 1:
		m_1 = ''
	if m_2 == 1:
		m_2 = ''
	if m_3 == 1:
		m_3 = ''
	
	title = r'$\vec{k} =$'
	
	if m_1 != 0:
		if len(title) != 11:
			title += r'$ +$'
		title += r'$ {} \vec{{b}}_1$'.format(m_1)
		
	if m_2 != 0:
		if len(title) != 11:
			title += r'$ +$'
		title += r'$ {} \vec{{b}}_2$'.format(m_2)
		
	if m_3 != 0:
		if len(title) != 11:
			title += r'$ +$'
		title += r'$ {} \vec{{b}}_3$'.format(m_3)
		
	elif m_1 == 0 and m_2 == 0 and m_3 == 0:
		title += r'$ \vec{{0}}$'
	
	ax.set_title(title, fontsize=20, y=1.08)

slid_amp.on_changed(update_phonon)
slid_m_1.on_changed(update_phonon)
slid_m_2.on_changed(update_phonon)
slid_m_3.on_changed(update_phonon)

# Freeze button
ax_freeze = plt.axes([0.8, 0.01, 0.1, 0.04])
button_freeze = Button(ax_freeze, 'Freeze', color=axcolor_1, hovercolor='0.975')

def freeze(event):

	slid_amp.reset()
	slid_m_1.reset()
	slid_m_2.reset()
	slid_m_3.reset()

	set_lattice_offsets(positions)
	excite_phonon(positions, radio_branch, [0, 0, 0], 0)
	
	ax.set_title('Cubic Lattice Phonons', fontsize=20, y=1.08)
    
button_freeze.on_clicked(freeze)

# Branch selector
ax_radio = plt.axes([0.025, 0.55, 0.2, 0.18], facecolor=axcolor_3)
button_radio = RadioButtons(ax_radio, ('Branch 1', 'Branch 2', 'Branch 3'), active=0)

def select_branch(label):

	global radio_branch
	radio_branch = int(label[-1])
	update_phonon(None)

button_radio.on_clicked(select_branch)

# Background toggle
ax_bgd = plt.axes([0.025, 0.40, 0.2, 0.12], facecolor=axcolor_3)
button_bgd = RadioButtons(ax_bgd, ('Frame On', 'Frame Off'), active=0)

def toggle_bgd(label):

	if label == 'Frame On':
		ax.axis('on')
	elif label == 'Frame Off':
		ax.axis('off')
    
button_bgd.on_clicked(toggle_bgd)

'''
Animate!
'''

ani = FuncAnimation(fig, update, tot_frames, interval=1000/fps, blit=False)
plt.show()

'''
TODO LIST:
	
	- SCRAP THE LCM THING
	- Just resize the position array everytime you wanna compute the mode
	
	- Precompute and store:
		- Original lattice offsets
	- Remove a, b, c (just set all = 1)
	- Make everything faster with numba (hop3fully)
	- Build bcc and fcc dynamical matrices
	- Consider diatomic systems
		- Diamond
		- Graphite
		- Zincblende

'''

'''
PHONON EQUATIONS DERIVATION:

For an orthogonal lattice, the dynamical matrix is given by:

\vec{D}( \vec{k} = k_x \vec{x} + k_y \vec{y} + k_z \vec{z} ) = 

\frac{4f}{m} * np.diagflat([np.square(np.sin(k_x * a /2)), \
							np.square(np.sin(k_y * b /2)), \
							np.square(np.sin(k_z * c /2))])
							
The frequencies are trivial:

\omega_1 = 2 * np.sqrt(f / m) * np.abs(np.sin(k_x * a /2))
\omega_2 = 2 * np.sqrt(f / m) * np.abs(np.sin(k_y * b /2))
\omega_3 = 2 * np.sqrt(f / m) * np.abs(np.sin(k_z * c /2))

The polarizations are:

\vec{v}_1 = \vec{x}
\vec{v}_1 = \vec{y}
\vec{v}_1 = \vec{z}

Finally, the atomic displacements are given by:

\vec{u}_i = \vec{v}_i * \exp{\vec{k} \cdot \vec{R}_n - \omega_i * t}

Here, \vec{k} as before and \vec{R}_n = n_1 \vec{a}_1 + n_2 \vec{a}_2 + n_3 \vec{a}_3
Thus, the atomic displacements, in terms of the real part, are:

\vec{u}_i = \vec{v}_i * \cos{k_x * n_1 * a + k_y * n_2 * b + k_z * n_3 * c + \tilde{t}}

We don't care about the specifics of the time, all we need to know is that
\tilde{t} varies from 0 to 2Pi to complete one cycle

Finally, we musn't forget that only certain vectors \vec{k} are permitted due to 
finite lattice sizes. We write:

\vec{k} = m_1 / N_1 * \vec{b}_1 + m_2 / N_2 *\vec{b}_2 + m_3 / N_3 * \vec{b}_3

Here, the N_i are the numbers of unit cells in each primitive direction.
The m_i integers take on the values m_i = 0, 1, ... , N_i - 1
Finally, using this form of \vec{k}, the atomic displacements can be written

\vec{u}_i = \vec{v}_i * \cos{ m_1 * n_1 / N_1 * 2Pi + \
							  m_2 * n_2 / N_2 * 2Pi + \
							  m_3 * n_3 / N_3 * 2Pi + \
							  - \tilde{t} }
							  
Again: 
- N_i fixed by crystal upon generation
- n_i determined by atomic location within crystal, n_i = 0, 1, ..., N_i - 1
- m_i chosen by specifying \vec{k}, m_i = 0, 1, ..., N_i - 1

We can save ourselves some multiplication effort if we finally define a final time
parameter \tau \in [0,1):

\vec{u}_i = \vec{v}_i * \cos{2Pi * (m_1 * n_1 / N_1 + \
									m_2 * n_2 / N_2 + \
									m_3 * n_3 / N_3 + \
									\tau)}
									
Final concern: the substitution \tau = \omega_i * t / 2Pi only makes sense in the
context of a single phonon. If we consider a superposition of different phonons of
DIFFERENT \omega_i's, then this does not work anymore, and then a larger period
must be chosen, which is a super cycle of each of the individual phonon cycles.

To deal with this, let us consider the frequency:
\omega_1 = 2 * np.sqrt(f / m) * np.abs(np.sin(k_x * a /2))

Immediately, scale the frequency in units of np.sqrt(f / m)

\tilde \omega_1 = 2 np.abs(np.sin(k_x * a /2))

What values of k_x can occur? Since k_x * \vec{x} = m_1 / N_1 * \vec{b}_1, and we know 
\vec{b}_1 = 2Pi / a \vec{x}, then

k_x = m_1 / N_1 * 2Pi / a

Thus

\tilde \omega_1 = 2 np.abs(np.sin( m_1 / N_1 * Pi ))

This is consistent with m_1 = 0, 1, ..., N_1 - 1

To highlight the different frequencies of the different phonons, we need to make the
length of frames required for one oscillation to be proportional to the frequency. 
Further, we need the position array to be of length that can be divided by the length
of frames required for each possible phonon. Since there are N_1 - 1 different 
frequencies (ignoring \omega=0) we need to find the lowest common multiple of these

num_frames_one_phonon = int( fps / omega * factor)

If we want the fastest oscillating phonon to have a period of T_fastest seconds, then

num_frames_one_phonon = int (fps / omega * omega_fastest * T_fastest)

To determine the length of the position array, we need

LCM ( num_frames( m_1=1 ), num_frames( m_1=2 ), ..., num_frames( m_1=N_1-1 ) )

where each num_frames is a function of omega, which are functions of m_1.
However, since N_1 need not equal N_2 need not equal N_3, we also LCM for all
three dimensions

How to compute this?

m_1 = np.arange(N_1 - 1) + 1
omega_1 = 2 * np.sin( m_1 / N_1 * Pi )
num_frames_1 = int( fps / omega_1 * omega_fastest * T_fastest )
god_frames_1 = np.lcm.reduce(num_frames_1)
...
god_frames = np.lcm.reduce(god_frames_1, god_frames_2, god_frames_3)
positions = np.empty(( 3 * num_atoms, god_frames ))

How does this affect the phonons?
The total number of frames required for one oscillation is

tot_frames = int( fps / omega_i * omega_max * T_max )

where omega_i \in [0, omega_max]

In the displacement vector, there is a cos(... - omega * t)

We focus on the omega * t term. To express t in terms of the given frame, let

t = \chi * current frame

s.t.

\omega_i * tot_frames * \chi = 2Pi

\chi = 2Pi / omega_i / tot_frames = 2Pi / fps / omega_max / T_max

-> \omega t = omega_i * (2Pi / fps / omega_max / T_max) * current_frame
'''
