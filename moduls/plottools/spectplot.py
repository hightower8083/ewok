#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from pylab import yscale, savefig, plot,xlim,ylim, clf, xlabel, ylabel
import numpy as np
from fnmatch import fnmatch
from os import listdir,system

system('mkdir -p img')

names = np.sort(listdir('./phase_data/'))

times_all = []

ext_x = (1,20)
#ext_y= (1.0e-2,10)
data_types = {'x':0,'y':1,'vx':2,'vy':3,'gamma':4}


for name in names:
        if fnmatch(name, '*000'):
                tim=name.split('_')[-2]
                times_all.append(tim)
print times_all, 'snaps are detected'

data_type = raw_input('input value to plot:')

times_in = raw_input('input time step:')

if times_in == '':
	times = times_all
else:
	times = str(times_in)
	while len(times)<4:
		times = '0'+times
	times=np.array([times])
	print times

for time in times:
	print 'doing '+time
	names_time = []
	for name in names:
		if fnmatch(name, 'phs_dat_'+time+'_*'):
			names_time.append('./phase_data/'+name)

	data = np.empty(0)
	for  name in names_time:
		aa = np.load(name)
		data = np.append(data, aa[data_types[data_type],:])

	hi, bins = np.histogram(data,bins=200, normed=True, range=ext_x)
	plot( 0.5*(bins[1:]+bins[:-1]), hi )
	yscale('log')
	xlim(ext_x)
#	ylim(ext_y)
	xlabel(r'$\gamma$', fontsize=18)
	ylabel(r'EDF [normed]', fontsize=18)
	savefig('./img/spec'+time+'.png')
	clf()


