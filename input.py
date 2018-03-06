import numpy as np
from read_opts import *

###################################
#### options evaluation  #########

courant = 0.98

dx = 1./steps
dy = dx
dt = dx*courant

neutral =  bool(neutral)
e_static = bool(e_static)

diag_phase = bool(diag_phs)
diag_field = bool(diag_fld)
diag_dens = bool(diag_dens)

profiles = ['cosinus','flat','gauss']
seed_profile = profiles[int(seed_profile)]
Nstat = int(num_part)

beta = np.sqrt(1-gamma_ext**-2)

l_x, l_y = beam_right-origin_x, beam_top-origin_y
bins_x, bins_y = np.array( np.rint(  [(-left+right)/dx , (top-bottom)/dy] ) ,dtype = 'int')

bins_y = int(2*np.ceil( ((top-bottom)/dy+1.)/2.)-1) #correction to keep bins_y even
top = bottom + bins_y*dy # correction to keep bins_y even

N_part_full = int(np.rint( Nstat*l_x*l_y/(dx*dy)  ))

N_max = int( np.rint(  Nstat*(right-left)*(top-bottom)/(dx*dy*num_node) ))

num_steps = int(np.ceil(simulation_time/dt))
injection_time = int(np.ceil(injection_time/dt))
snapshots_interval = int(np.ceil(snapshots_interval/dt))

seed_time = int(seed_time//dt)
seed_delay = int(seed_delay//dt)
seed_duration = seed_time + seed_delay

#seed_duration = int(np.ceil(seed_duration/dt))

###################################
####  manually added options ######

e_charge = -1.0
M_p = 1836.1253211

density_profile = 'homogen'  # known are 'homogen', 'ramp_lin', 'exp1', 'circle'
prof_length = 0.9 # sclae length of density profile (ramp, rms width, radius)

EDFx= 'thermal'  # phase distribution may either 'thermal' or 'waterbag'
EDFy= 'thermal'

stepwise_step = 3      # spatial step in number of cells used for field and dens stepwise slice output

#_________ particle_tracking ______________

tracking = False
marked_veloc = np.array([[0.01,0.01]])
marked_coord = np.array([[8, 0.]])

markID = 3 + N_max*int(marked_coord[0,0]/(-left+right)*num_node)

#________signal asymetry (flat case)_______

asym = False

#___________boundary options_______________

x_period = True # periodicity for e.m. field in x-direction
y_period = True # periodicity for e.m. field in y-direction
pillow_size = 20  # length of absorbsion lines at transvers boundaries, in cell numbers (only when y_period=False)

x_particle_absorb = False
y_particle_absorb = False

#_______oblique_signal_injection___________

seed_angle = 0.0

#_________pump_modifications_______________
bonus = False
a_bonus = 0.
w_bonus = 0.

lattice_modif = False # disbalanced pump arms

pump_unstable = False  # oscillation of pump amplitude
delta_ampl = 0.45
delta_freq = 1.e-2 # in units of omega_0
##################################
