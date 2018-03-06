#! /usr/bin/env python
import numpy as np
from mpi4py import MPI
from sys import path
path.append('./moduls/')
from ewoks import *
from input import *

comm = MPI.COMM_WORLD
ewok = Ewoks(comm)
time_steps = np.arange(0,num_steps+1)

for time in time_steps:
    ewok.toy.diagnostics_stepwise(ewok,time)
    ewok.toy.diagnostics_full(ewok,time)

    ewok.potential_calcul(time)

    ewok.motion()

    if e_static and neutral:
        ewok.motion_proton()

ewok.toy.bye_msg(ewok)
