import numpy as np
#######################################
#### reading options from text file ###
a0, \
injection_time,\
a_seed,\
seed_time,\
seed_delay,\
w_seed,\
seed_profile,\
n_e,\
gamma_ext,\
v_x,\
v_y,\
beam_top,\
origin_x,\
beam_right,\
origin_y,\
num_part,\
simulation_time,\
steps,\
top,\
left,\
right,\
bottom,\
snapshots_interval,\
diag_phs,\
diag_fld,\
diag_dens,\
num_node,\
seed_width,\
e_static,\
neutral,\
out_folder = np.loadtxt('./input.txt',dtype='str')

notification = False
user_email = ["ewok@dark_side.ndr"]

a0,injection_time, a_seed,seed_time,seed_delay,w_seed,\
seed_profile,n_e,gamma_ext,v_x,v_y,beam_top,\
origin_x,beam_right,origin_y,num_part,simulation_time,steps,\
top,left,right,bottom,snapshots_interval,diag_phs,\
diag_fld,diag_dens,num_node,seed_width,e_static,neutral = np.array([\
a0,injection_time,a_seed,seed_time,seed_delay,w_seed,\
seed_profile,n_e,gamma_ext,v_x,v_y,beam_top,\
origin_x,beam_right,origin_y,num_part,simulation_time,steps,\
top,left,right,bottom,snapshots_interval,diag_phs,\
diag_fld,diag_dens,num_node,seed_width,e_static,neutral],dtype='d')
