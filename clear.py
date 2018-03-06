#! /usr/bin/env python
from os import system
from input import out_folder

system('rm -rf ' + out_folder)
system ('mkdir '+out_folder)
system ('mkdir '+out_folder+'phase_data/')
system ('mkdir '+out_folder+'phase_ion/')
system ('mkdir '+out_folder+'fld_data/')
system ('mkdir '+out_folder+'dens_data/')
system ('mkdir '+out_folder+'dens_ion/')
system ('mkdir '+out_folder+'stat_potent/')
system ('mkdir '+out_folder+'stepwise_data/')
system('cp -f input.py input.txt read_opts.py ./moduls/plottools/plotter.py '+out_folder)
