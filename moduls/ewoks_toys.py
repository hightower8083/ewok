import numpy as np
from mpi4py import MPI
import time as timmy
from scipy.fftpack import fft, ifft
import os
from input import *

text_out1 = '**Simulation is launched on %s  nodes. \nYou have %s electrons spread on %s x %s grid. \nInitializing ewoks objects.\n' 
text_out2 = "**Cleaning the output folders\n"
text_out3 = "**All initial data is calculated. Starting the main simulation loop..\n"
text_out4 = "** %s percent is done. output files are written.\n"
text_out5 = "**Im done in %s minutes. Bye !\n"

exp_src = np.vectorize( lambda x, y, t_step : np.exp(-2*np.pi*1j*((t_step-seed_delay)*dt-x)))

class toys:
    def __init__(self, ewok):
        '''gives Ewok his toys for signal and pump creation, diagnostics and messaging'''
        self.x_ij, self.y_ij = ewok.x_ij, ewok.y_ij

        # preparing the optical lattice potential
        if lattice_modif:
            self.stat_pump_potent = -a0*(1./2.)*1j*np.exp(-2j*np.pi*((beta+1)*self.x_ij + self.y_ij/gamma_ext)) \
              + a0*(3./2.)*1j*np.exp(-2j*np.pi*((beta+1)*self.x_ij - self.y_ij/gamma_ext))
        else:
            self.stat_pump_potent = 2*a0*np.sin(2*np.pi*self.y_ij/gamma_ext)\
                * np.exp(-1j*(beta+1)*self.x_ij*2*np.pi)

        # preparing the seed
        self.get_seed_profile()
        self.seed_rotation()

        # preparing the working folders, writing message and starting the timer
        if ewok.rank == 0:
            os.system('./clear.py')
            text_out = open(out_folder+'text_out.txt',mode='a')
            text_out.writelines(text_out2)
            text_out.writelines(text_out1 % (ewok.size, N_part_full,bins_x,bins_y))
            text_out.close()
            self.Tstart = timmy.time()

        # defining the output folders names
        self.stepwise_folder = out_folder+'stepwise_data/'
        self.phase_folder = out_folder+'phase_data/'
        self.ion_phase_folder = out_folder+'phase_ion/'
        self.dens_folder = out_folder+'dens_data/'
        self.ion_folder = out_folder+'dens_ion/'
        self.fld_folder = out_folder+'fld_data/'
        self.stat_potent_folder = out_folder+'stat_potent/'

    def pump_func(self,t):
        ''' temporal modulations of lattice potential '''
        pump =  self.stat_pump_potent

        if t<=injection_time:
            pump = pump*float(t)/float(injection_time)

        if bonus:
            pump = pump+self.bonus_pump(t)

        if pump_unstable:
            pump = pump*(1.0 + delta_ampl*np.sin(delta_freq*2*np.pi*t*dt))

        return pump

    def bonus_pump(self, time):
        ''' additional wave in the lattice (not used)'''
        if time<injection_time:
            return a_bonus*(time*injection_time**-1)*np.exp(-1j*(w_bonus-1)*(time+self.x_ij)*2*np.pi)
        else:
            return a_bonus*np.exp(-2j*np.pi*(w_bonus-1)*(time+self.x_ij))

    def seed(self,time):
        ''' constructs the Fourier image of the signal with defined frequency and profile'''

        if seed_profile == 'gauss' or seed_profile == 'cosinus':
            temporal = np.sin(np.pi*(time-seed_delay)/(seed_duration-seed_delay))
        elif seed_profile == 'flat':
            temporal = 1.

        seed =  temporal * a_seed * self.profile * self.rotator \
            * np.exp(-2*np.pi*1j*(w_seed - 1)*(time-seed_delay)*dt - np.pi*1j/2.)

        return fft(seed)

    def seed_rotation(self):
        ''' constructs the rotation multiplier for the signal (doesn't turn the phase front)'''
        if seed_angle == 0.0:
            Ks_perp = 0.0
        else:
            Ks_perp = (1+seed_angle**-2.)**-0.5

        yy = np.r_[top:bottom:(bins_y+1)*1j]
        self.rotator =  np.exp(-2*np.pi*1j*Ks_perp*yy)

    def diagnostics_full(self,ewok,time):
        '''manager-script for main 2D diagnostics: density, phase coordinates, e.m. field, e.static potential'''
        if snapshots_interval!= 0 and np.mod(time,snapshots_interval)==0.0:
            self.data_output(ewok,time)
            if ewok.rank == 0:
                text_out = open(out_folder+'text_out.txt',mode='a')
                precent_done = str(np.int(np.rint(100*time/num_steps)))
                text_out.write(text_out4 % precent_done)
                text_out.close()

    def diagnostics_stepwise(self,ewok,time):
        ''' writing stepwise diagnosics: field and density dynamics at left boundary; probe particle track(trigs the method)'''
        if tracking:
            self.track_out(ewok)

        if ewok.rank == ewok.size-1:
            fld_tmp = ewok.potent_em[[-1],::stepwise_step]
            exp_tmp = exp_src(self.x_ij, self.y_ij,time)[[-1],::stepwise_step]
            file_fld = open(self.stepwise_folder+'field_oscil', mode='a')
            file_dens = open(self.stepwise_folder+'density_amplitude', mode='a')

            np.savetxt(file_fld, np.real(fld_tmp*exp_tmp))
            np.savetxt(file_dens, \
                np.real(ewok.dens2gamma[[ewok.bins[0]//2],::stepwise_step]))

            file_fld.close()
            file_dens.close()

    def data_output(self, ewok,time):
        ''' creates 2D data output (main output method)'''

        me = str(ewok.rank)
        while len(me) < 3: me = '0'+me

        tt = str(int(time*dt))
        while len(tt) < 4: tt = '0'+tt

        if diag_phase:
            phase_file = self.phase_folder+'phs_dat_'+me+'_'+tt
            np.save(phase_file, ewok.particles )
            if e_static and neutral:
                phase_file = self.ion_phase_folder + 'iphs_dat_'+me+'_'+tt
                np.save(phase_file,ewok.particles_p)

        if diag_dens:
            ewok.dens_calcul()
            data_pass = np.array(ewok.data, dtype='d')

#            ewok.dens2gamma_calcul()
#            data_pass = np.array(ewok.dens2gamma, dtype='d')

            if ewok.rank==0:
                data_tot = np.empty(( (ewok.bins[0]+1)*ewok.size, ewok.bins[1]+1 ),dtype='d')
            else:
                data_tot = None

            ewok.comm.Gather([data_pass,MPI.DOUBLE],[data_tot,MPI.DOUBLE], root=0)

            if ewok.rank == 0:
                data_tot = np.delete(data_tot,ewok.to_kick,0)
                dens_file = self.dens_folder+'dens_dat_'+tt
                np.save(dens_file, data_tot)

        if diag_field:
            potent_pass = np.array(np.real(ewok.potent_em*exp_src(ewok.x_ij, ewok.y_ij,time)),dtype='d')
            if ewok.rank == 0:
                data_tot = np.empty(( (ewok.bins[0]+1)*ewok.size,ewok.bins[1]+1 ),dtype='d')
            else:
                data_tot = None

            ewok.comm.Gather([potent_pass,MPI.DOUBLE],[data_tot,MPI.DOUBLE], root=0)

            if ewok.rank == 0:
                data_tot = np.delete(data_tot,ewok.to_kick,0)
                vect_potent_file = self.fld_folder+'fld_dat_'+tt
                np.save(vect_potent_file, data_tot)

        if e_static and diag_field and ewok.rank == 0:
            stat_potent_file = self.stat_potent_folder+'stat_potent_'+tt
            np.save(stat_potent_file,np.real(ewok.stat_potents))

        if e_static and neutral and diag_dens:
            data_pass = np.array(ewok.dens_p,dtype='d')

            if ewok.rank==0:
                data_tot = np.empty(( (ewok.bins[0]+1)*ewok.size,ewok.bins[1]+1 ),dtype='d')
            else:
                data_tot = None

            ewok.comm.Gather([data_pass,MPI.DOUBLE],[data_tot,MPI.DOUBLE], root=0)

            if ewok.rank == 0:
                data_tot = np.delete(data_tot,ewok.to_kick,0)
                dens_file = self.ion_folder+'idens_dat_'+tt
                np.save(dens_file, data_tot)

    def track_out(self,ewok):
        ''' writing probe particle track (actual method)'''
        indx = np.where(ewok.particles[:,-1] ==markID)[0]
        if indx.shape[0]!=0:
            track = open(self.stepwise_folder+'particle_track',mode='a')
            np.savetxt(track, ewok.particles[indx,:])
            track.close()

    def get_seed_profile(self):
        '''creates the seed transverse spatial profile'''

        seed_axis = np.arange(bins_y+1)
        relat_width = (seed_axis*dy-0.5*(top-bottom))/seed_width
        profile_cutout = np.array(np.abs(relat_width)<=1.,dtype=int)

        if seed_profile == 'flat':
            if asym:
                self.profile= np.sin(3.*np.pi*relat_width*profile_cutout)
            else:
                self.profile = profile_cutout
        elif seed_profile == 'cosinus':
            self.profile = np.cos(np.pi/2.*relat_width)*profile_cutout
        elif seed_profile == 'gauss':
            self.profile = np.exp(-relat_width**2)
        else:
            print( 'seed profile is not understood')


    def start_msg(self,ewok):
        '''text output message before mainloop start'''
        if ewok.rank == 0:
            text_out = open(out_folder+'text_out.txt',mode='a')
            text_out.writelines(text_out3)
            text_out.close()

    def bye_msg(self,ewok):
        '''text output message after mainloop with timing and email notification '''
        ewok.Npart_absorbed = ewok.comm.reduce(ewok.Npart_absorbed, root=0)
        if ewok.rank == 0:
            text_out = open(out_folder+'text_out.txt',mode='a')
            text_out.write( text_out5 % str((timmy.time()-self.Tstart)/60))
            text_out.write ( str(ewok.Npart_absorbed)+' particles are absorbed at transverse boundaries')
            text_out.close()
