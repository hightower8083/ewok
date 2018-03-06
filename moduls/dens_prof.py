import numpy as np
from input import *

#####  profile functions #####

def ramp_lin(x):
    if x<=prof_length:
        return (x/prof_length)
    else:
        return 1.0

ramp_lin = np.vectorize(ramp_lin)

ramp_exp = lambda x: np.exp( -np.abs( 2*(x-left - prof_length)/prof_length   )**2.5 )

def circ(x,y):
    if np.sqrt((x- 0.5*(right-left))**2+y**2)<=prof_length**2:
        return 1.0
    else:
        return 0.

circ = np.vectorize(circ)

#########################

def profiled(ewok):
    if density_profile == 'homogen':
        ewok.particles[:,-2] = np.ones(ewok.particles.shape[0], dtype='d')
    elif density_profile == 'ramp_lin':
        ewok.particles[:,-2] = ramp_lin(ewok.particles[:,0])
    elif density_profile == 'exp1':
        ewok.particles[:,-2] = ramp_exp(ewok.particles[:,0])
    elif density_profile == 'circle':
        if ewok.particles.shape[0]!=0:
            ewok.particles[:,-2] = circ(ewok.particles[:,0], ewok.particles[:,1])
    ewok.particles[:,-2] = ewok.max_wght * ewok.particles[:,-2]

def particle_remover(ewok):
    if (ewok.particles[:,-2] == 0.0).sum() !=0:
        particles2keep_indx = np.nonzero( ewok.particles[:,-2] != 0.0 )[0]
        ewok.particles = ewok.particles[particles2keep_indx,:]
