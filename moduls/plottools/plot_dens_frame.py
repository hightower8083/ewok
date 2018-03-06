#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from pylab import imshow,colorbar,savefig,clf, xlim, ylim
from fnmatch import fnmatch
from input import *
from os import listdir, system

system('mkdir -p img')

names = np.sort(listdir('./dens_data/'))

times_all = []

for name in names:
	tim=name.split('_')[-1]
	times_all.append(tim)

print times_all, 'snaps are detected'


xmin, xmax,ymin,ymax, vmin,vmax = np.array(raw_input('input limits (xmin, xmax,ymin,ymax,vmin,vmax):').split(','),dtype='float')

times_in = raw_input('input time step:')


if times_in == '':
	times = times_all
	print 'writing snaps',times
else:
	while len(times_in)<4:
		times_in = '0'+times_in
	times=np.array([times_in])

for time in times:
	if int(time)>=xmax:
		data = np.load('./dens_data/dens_dat_'+time)
		ext=np.array([left,right,bottom,top])
		if vmin==vmax:
			vmin=vmax=None
		imshow(data.transpose(),origin='lower',aspect='auto',extent=ext, cmap='spectral',interpolation='spline36',vmin=vmin,vmax=vmax);colorbar()
		xlim(int(time)-xmax, int(time))
		if ymin!=ymax:
			ylim(ymin,ymax)
		savefig('./img/dens_'+time+'.png')
		clf()
		print 'file is written'
