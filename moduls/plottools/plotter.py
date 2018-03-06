#! /usr/bin/env python

from matplotlib import use
import readline
import numpy as np
from fnmatch import fnmatch
from input import *
from os import listdir, system
from sys import argv


system('mkdir -p img') 

out_fold = {'phs':'./phase_data/', 'fld': './fld_data/', 'dns': './dens_data/', 'estat': './stat_potent/',\
            'iphs':'./phase_ion/', 'idns': './dens_ion/','chrg':'./dens_data/'}
vals = {'x':0,'y':1,'vx':2,'vy':3,'gamma':4,'wgt':5}
labels = {'x':r'$x/\lambda_0$','y':r'$y/\lambda_0$','vx':r'$v_x/c$','vy':r'$v_y/c$','gamma':r'$\gamma_e$'}
file_mask =  {'phs':'phs_dat_', 'fld': 'fld_dat_', 'dns': 'dens_dat_', 'estat': 'stat_potent_', \
              'idns': 'idens_dat_', 'iphs':'iphs_dat_','chrg':'chrg_dat_'}
img_format='png'

class plotter:
	def __init__(self):
		if len(argv)>2:
			self.corient = argv[1]
			self.formt='%g'
		else:	
			self.corient = 'vertical'
			self.formt=None
		self.cmap = 'spectral'
		reg = argv[-1]
		if reg == 'i':
			self.interactive = True
			from pylab import ion
			ion()
		else:
			self.interactive = False
			use('Agg')
		from pylab import figure
		self.fig = figure()
		while True:
			dat_reg = raw_input('choose type of data to plot:').split(' ')			
			if len(dat_reg)>1 and dat_reg[-1]!='s':
				self.cmap = dat_reg[-1]
#			else:   self.cmap = 'spectral'
			self.val_name  = dat_reg[0]
			if self.val_name=='': return
			if self.val_name in ('fld', 'dns','estat', 'idns','chrg'):
				if 's' in dat_reg:
					doer = self.grid_spec_plotter
				else:
					doer = self.grid_plotter
			elif self.val_name in ('phs', 'iphs'):
					doer = self.phase_plotter
			elif self.val_name == 'help':
				print 'AAA!! Please somebody help!'
				doer = lambda : 1
			self.names = np.sort(listdir( out_fold[self.val_name] ))
			self.time_detect()
			self.time_prompt()
			doer()

	def phase_data_read(self,time):
                time = str(time)
                while len(time)<4:
	        	time = '0'+time
                names_time = []
                for name in self.names:
	        	if fnmatch(name, file_mask[self.val_name]+'*_'+time):
                  		names_time.append(out_fold[self.val_name]+name)
        	data1 = np.empty(0)
	        data2 = np.empty(0)
		wght = np.empty(0)
                for name in names_time:
	        	dat = np.load(name,mmap_mode='r')
        	        data1 = np.append( data1, dat[:,vals[self.type1]] ) 
                	data2 = np.append( data2, dat[:,vals[self.type2]] )
			wght = np.append( wght, dat[:,vals['wgt']] )
               	self.data,self.ax1,self.ax2 = np.histogram2d(data1,data2,bins=[self.bins_x,self.bins_y],\
		     weights = wght, range=[self.ext1,self.ext2])

	def phase_plotter(self):
		from pylab import imshow, axis,show,colorbar,savefig,clf,ion,figure,ioff,xlabel,ylabel
		self.type1, self.type2  = raw_input('input axes:').split(',')
		plot_domain = raw_input('input limits (xmin, xmax,ymin,ymax,vmin,vmax):').split(',')
                if plot_domain[0]=='':
                        xmin = xmax = ymin = ymax = vmin = vmax = 0
                else: xmin, xmax,ymin,ymax, vmin,vmax = np.array(plot_domain,dtype='float')
		if xmin == xmax:
			self.ext1 = self.ext_def(self.type1)
                else:	self.ext1 = [xmin, xmax]
                if ymin == ymax:
			self.ext2 = self.ext_def(self.type2)	
                else:	self.ext2 = [ymin,ymax]
		if vmin==vmax:
                	vmin=vmax=None
		bins_raw = raw_input('input bins:').split(',')
		if bins_raw[0]=='': 
			self.bins_x, self.bins_y = bins_x, bins_y
		else:	self.bins_x, self.bins_y = np.array(bins_raw, dtype ='i')
		
	        for time in self.times:
			self.phase_data_read(time)
                	self.fig.clf()
			imshow(self.data.transpose(),extent=(self.ax1[0],self.ax1[-1],self.ax2[0],self.ax2[-1]), \
			   aspect='auto',origin='lower',cmap=self.cmap,interpolation='spline36',vmin=vmin, vmax=vmax)
			colorbar(orientation=self.corient,format=self.formt)
			xlabel(labels[self.type1],fontsize=18)
			ylabel(labels[self.type2],fontsize=18)
                	axis(self.ext1 + self.ext2)
			if self.interactive:
				self.fig.show()
			savefig('./img/'+file_mask[self.val_name]+str(self.type1)+str(self.type2)+'_'+str(time)+'.'+img_format)
			print 'file is written'

	def grid_plotter(self):
		from pylab import imshow,show,colorbar,savefig,xlim,ylim,xlabel,ylabel
		plot_domain = raw_input('input limits (xmin, xmax,ymin,ymax,vmin,vmax):').split(',')
		if plot_domain[0]=='':
			xmin = xmax = ymin = ymax = vmin = vmax = 0
		else: xmin, xmax,ymin,ymax, vmin,vmax = np.array(plot_domain,dtype='float')
		for time in self.times:
			while len(time)<4:
                        	time = '0'+time
			if self.val_name == 'chrg':
				self.data = np.load(out_fold['idns']+file_mask['idns']+time)-\
				            np.load(out_fold['dns']+file_mask['dns']+time)
			else:
	        		self.data = np.load(out_fold[self.val_name]+file_mask[self.val_name]+time)
        		ext=np.array([left,right,bottom,top])
        		if vmin==vmax:
		                vmin=vmax=None
			self.fig.clf()	
	        	imshow(self.data.transpose(),origin='lower',aspect='auto',extent=ext, \
				cmap=self.cmap,interpolation='spline36',vmin=vmin, vmax=vmax)
			colorbar(orientation=self.corient,format=self.formt)
			xlabel(labels['x'],fontsize=18)
			ylabel(labels['y'],fontsize=18)
	        	if xmin != xmax:
	                	xlim(xmin, xmax)
        		if ymin!=ymax:
		                ylim(ymin,ymax)
			if self.interactive:
                                self.fig.show()
	        	savefig('./img/'+file_mask[self.val_name]+time+'.'+img_format)
        		print 'file is written'

	def grid_spec_plotter(self):
                from pylab import imshow,show,colorbar,savefig,xlim,ylim

                plot_domain = raw_input('input limits (xmin, xmax,ymin,ymax,vmin,vmax):').split(',')
                if plot_domain[0]=='':
                        xmin = xmax = ymin = ymax = vmin = vmax = 0
                else: xmin, xmax,ymin,ymax, vmin,vmax = np.array(plot_domain,dtype='float')

                for time in self.times:
                        while len(time)<4:
                                time = '0'+time
			data = np.load(out_fold[self.val_name]+file_mask[self.val_name]+time)
			data_fft = np.abs(np.fft.fft2(data[:2*(data.shape[0]//2),:2*(data.shape[0]//2)]))
			data_fft = 4*data_fft[:data_fft.shape[0]//2,:data_fft.shape[1]//2]/(data_fft.shape[0]*data_fft.shape[1])
			ext = [0, np.fft.fftfreq(data.shape[0],dx).max(),\
			       0, np.fft.fftfreq(data.shape[1],dy).max()]
			if vmin==vmax:
                                vmin=vmax=None
			self.fig.clf()
			imshow(np.abs(data_fft).transpose(),aspect='auto', \
                		cmap=self.cmap,origin='lower',extent=ext,\
                		interpolation='spline36',vmin = vmin, vmax=vmax)
			colorbar(orientation=self.corient,format=self.formt)
                        if xmin != xmax:
                                xlim(xmin, xmax)
                        if ymin!=ymax:
                                ylim(ymin,ymax)
                        if self.interactive:
                                self.fig.show()
                        savefig('./img/'+file_mask[self.val_name]+'spec_'+time+'.'+img_format)
                        print 'file is written'

	def ext_def(self,type):
		if type == 'x':
			return [left,right]
		elif type == 'y':
			return [bottom,top]
		elif type == 'vx':
			return [-2*v_x,2*v_x]
		elif type == 'vy':
			return [-2*v_y,2*v_y]

	def time_detect(self):
		self.times_all = []
		for name in self.names:
			if self.val_name in ('phs', 'iphs'):
				if fnmatch(name, '*_000_*'):
	        			tim=name.split('_')[-1]
        				self.times_all.append(tim)
			else:	
				tim=name.split('_')[-1]
				self.times_all.append(tim)
		print self.times_all, 'snaps are detected'
		
	def time_prompt(self):
        	times_in = raw_input('input time step:').split(',')
        	if times_in[0] == '':
	                self.times = self.times_all
        	else:
	                self.times = times_in

## construct/run the plotter
plotter()
#######################
