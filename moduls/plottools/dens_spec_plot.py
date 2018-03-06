from pylab import *
from fnmatch import fnmatch
from os import listdir

bottom,top = -2.5, 2.5

names = listdir('./')

dens_names = []
for name in names:
	if fnmatch(name, 'density_*'):
		dens_names.append(name)
dens_names.sort()

fld_names = []
for name in names:
        if fnmatch(name, 'field_*'):
                fld_names.append(name)
fld_names.sort()

for name in dens_names:
	a3 = np.load(name)
	dens_spect = np.abs(np.fft.fft(a3,axis=0))[:a3.shape[0]//2,:]/a3.shape[0]
	imshow(dens_spect.transpose(),aspect='auto',cmap='jet',origin='lower',extent=[0,12.5,bottom,top],interpolation='spline36',vmax=3e-6)
	xlim(1e-3,0.033)
	xlabel(r'$\omega_e/\omega_0$',fontsize=18)
	ylabel(r'$x/\lambda_0$',fontsize=18)
	savefig(name+'_img.png')
	clf()
	print 'density spectrum written'

for name in fld_names:
        a2 = np.load(name)
	fld_spect = np.abs(np.fft.fft(a2,axis=0))[:a2.shape[0]//2,:]/a2.shape[0]
	imshow(fld_spect.transpose(),aspect='auto',cmap='jet',origin='lower',extent=[0,12.5,bottom,top],interpolation='spline36',vmax=1.5e-5)
	xlim(0.97,1.03)
        xlabel(r'$\omega_s/\omega_0$',fontsize=18)
        ylabel(r'$x/\lambda_0$',fontsize=18)
        savefig(name+'_img.png')
        clf()
        print 'field spectrum written'
