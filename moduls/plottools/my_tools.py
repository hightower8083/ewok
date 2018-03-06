from pylab import *
from os import listdir
from fnmatch import fnmatch
#from scipy import *

vals = {'x':0,'y':1,'vx':2,'vy':3,'gamma':4,'wgt':5}


class tools():


	def averaging(self, d, dx = 0.04, avrg_step = 4):
		d = 2*d**2
                tim = dx*np.arange(d.shape[0])
                dstep = avrg_step*int(1./dx)
                a_avrg = np.zeros(0)
                coord_avrg = np.zeros(0)
                for i in np.arange(0,d.shape[0]-dstep, dstep):
                        a_avrg = np.hstack(( a_avrg, np.sqrt(d[i:i+dstep].mean()) ))
                        coord_avrg = np.hstack(( coord_avrg, tim[i:i+dstep].mean() ))
                return coord_avrg, a_avrg

        def signal_averaging(self, folder, dt=0.04, line=None, avrg_step = 4):
		d = np.load(folder+'stepwise_data/field_oscil',mmap_mode='r')
                if line==None:
                        line = d.shape[1]//2+1
#                d = 2*np.abs(d[:,line])
		d = 2*d[:,line]**2
		tim = dt*np.arange(d.shape[0])
		dstep = avrg_step*int(1./dt)
		a_avrg = np.zeros(0)
		tim_avrg = np.zeros(0)
		for i in np.arange(0,d.shape[0]-dstep, dstep):
			a_avrg = np.hstack(( a_avrg, np.sqrt(d[i:i+dstep].mean()) ))
			tim_avrg = np.hstack(( tim_avrg, tim[i:i+dstep].mean() ))
		return tim_avrg, a_avrg


        def phase_data_read(self, filt=None):
		names = np.sort(listdir('./phase_data/'))
		type1, type2  = raw_input('input axes:').split(',')
                time = raw_input('input time step:')
                while len(time)<4:
                        time = '0'+time
                names_time = []
                for name in names:
                        if fnmatch(name, 'phs_dat_*_'+time):
                                names_time.append('./phase_data/'+name)
                self.data1 = np.empty(0)
                self.data2 = np.empty(0)
                self.wght = np.empty(0)
                for name in names_time:
                        dat = np.load(name,mmap_mode='r')
			if filt != None:
#				nrg = 0.5*(dat[vals['vx'],:]**2+dat[vals['vy'],:]**2)
#				filt_indx = np.nonzero( (nrg>filt[0])*(nrg<filt[1]))[0]
				vx = dat[vals['vx'],:]
				filt_indx = np.nonzero( (vx>filt[0])*(vx<filt[1]))[0]
				self.data1 = np.append( self.data1, dat[vals[type1], filt_indx] )
				self.data2 = np.append( self.data2, dat[vals[type2], filt_indx] )
				self.wght = np.append( self.wght, dat[vals['wgt'], filt_indx] )
			else:
                        	self.data1 = np.append( self.data1, dat[vals[type1],:] )
                        	self.data2 = np.append( self.data2, dat[vals[type2],:] )
                        	self.wght = np.append( self.wght, dat[vals['wgt'],:] )
			print '.',


	def growth_rate(self,d,start=2000,end = 25000, dt=0.04):
		xdata = dt*np.arange(start,end)
		ydata = log(abs(d[start:end]))
		coeffs = polyfit(xdata,ydata,1)
		return coeffs[0]/(2*pi)

        def growth_rate_fit(self,d,start=2000,end = 25000,  dt=0.04):
                xdata = dt*np.arange(start,end)
                ydata = log(abs(d[start:end]))
                coeffs = polyfit(xdata,ydata,1)
		yfit = polyval(coeffs, xdata)
                return xdata, yfit

	def growth_rate_avrgd(self,folder, line=None, start=200,end = 600, dt=0.04):
		d = np.load(folder+'stepwise_data/field_oscil',mmap_mode='r')
		if line==None:
			line = d.shape[1]//2+1
		d = d[:,line]
		start, end = np.array([start/dt, end/dt], dtype = 'int' )
		delta_window = (start - end)//4.
		step = delta_window//10.
		window = np.arange(-delta_window,delta_window,step)
		g = []
		for st in start + window:
			for en in end + window:
				g.append(tool.growth_rate(np.abs(d),start=st,end = en))
		g = np.array(g)
		return g.mean(), g.std()

	def all_detectors_read(self,index, nnodes = 8, address = '../XFEL_SCRATCH/xfel_', extension = '/etc/fld_'):
		d=[]
		index = str(index)
		nodes = np.arange(nnodes)
		while len(index)<2:
                	index = '0'+index
		for node in nodes:
			node = str(node)
			while len(node)<3:
                        	node = '0'+node
			d.append(loadtxt(address+index+extension+node+'_000'))
		return d

	def all_detectors_shaper(self,d):
		shape0 = np.shape(d)
		nnodes = shape0[0]
		detector_data = np.arange(nnodes)
		time_data = d[0][:,0]/2.5
		data = zeros((shape0[1],0))
		for i in detector_data:
			data = np.hstack((data,d[i][:,1:]))
		return data.transpose(), time_data, detector_data

	def arrays_read(self,index_array=[1,2,3,4,5,6,7,8],address = '../XFEL_SCRATCH/xfel_', extension = '/etc/fld_004_000'):
		"""
		reading modif-PICLS2D em field data from address folders with the names containig indexes index_array
		"""
		d=[]
		for i in index_array:
			indx_name = str(i)
			while len(indx_name)<2:
				indx_name = '0'+indx_name
			d.append(loadtxt(address+indx_name+extension))
		return d

tool = tools()


