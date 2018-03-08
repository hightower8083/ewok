import numpy as np
from mpi4py import MPI

from input import *
from scipy.fftpack import fft, ifft

import fort_ewok
from ewoks_toys import *
from density import *
from poissoner import *
from dens_prof import *

class Ewoks(density_rout):
    def __init__(self,comm):
        ''' creates the ewok object, which has: communicator, grid geometry data, empty arrays for potential,
		    parcticles phase coordinates, density data on grid, charge neutralization data and some utility data.
        '''
        if comm.rank == 0:
            print('initializing ewok')
        self.comm = comm	# passing MPI intracommunicator to Ewok

        self.get_geom()		# 2D grid initiate
        if comm.rank == 0:
            print('geom')

        self.wave_nums()	# define the coordinate axis for k_y

        self.particle_init()	# initiate the electrons with equal weights
        if comm.rank == 0:
            print('particle_init')

        profiled(self)		# assigning weights according to chosen density profile

        particle_remover(self)

        if comm.rank == 0:
            print('particle_remover')

        if e_static and neutral:
            self.proton_init()	# creating neutralizing protons
            if comm.rank == 0:
                print('proton_init')
        else:
            self.dens_p = 0.

        self.toy = toys(self)	# constructing TOYs object
        if comm.rank == 0:
            print('toys')

        self.potent_em  = np.zeros(self.bins+1,dtype='complex')		#
        self.potent_em_fft = np.zeros(self.bins+1,dtype='complex')	#
        self.potent_next = np.zeros(self.bins[1]+1,dtype='complex')	#
        self.pump_potent_stat = self.toy.stat_pump_potent   		# reassigning static lattice potential

        self.Npart_absorbed = 0

        if not y_period:
            self.tansverse_env()		# initiating horizontal absorption lines

        if self.rank == 0 and e_static:
            self.solver  = laplace_splu(self.bins_full)
            print('constructed solver with LU decomposition of Laplace operator')

        self.comm.barrier()

        self.leap_frog_motion()

        self.toy.start_msg(self)
        if self.rank == 0:
            print('simulation is started; for the progress see text_out')


    def get_geom(self):
        '''creates the local ewok grid geometry and data according to communicator values.'''
        # getting the number of calculation nodes
        self.size = self.comm.Get_size()

        # getting the self index of the node (self-identification)
        self.rank = self.comm.Get_rank()

        # domain decomposition(x-axis)
        self.bins = np.array([np.rint((right-left)/self.size/dx),bins_y],dtype='int')

        # adding cell-column, common with next-to-the-right node
        self.bins_full = np.array([bins_x, bins_y])+1

        # defining local vertical boundaries
        self.left_p = left + self.rank*(self.bins[0])*dx
        self.right_p = self.left_p + (self.bins[0])*dx

        # horizontal boundaries are common
        self.bottom_p, self.top_p = bottom, top

        # utility array for density arrays collecting
        self.to_kick = (self.bins[0]+1)*np.arange(1,self.size,1)

        # grid nodes coordinates
        self.x_ij, self.y_ij = np.mgrid[ self.left_p:self.right_p:(self.bins[0]+1)*1.0j,\
                                         self.bottom_p:self.top_p:(self.bins[1]+1)*1.0j ]

    def particle_init(self):
        ''' creates initial phase coordinates, IDs of the particles and tracked particle'''
        # defining maximum particle weight
        if Nstat != 0:
            self.max_wght = (n_e*dx*dy/Nstat)
        else:
            self.max_wght = 0.0

        # next goes initiation of electrons coordinates within the rectangular domain
        if self.left_p >= l_x+origin_x or self.right_p<=origin_x:
            self.particles =  np.zeros((0,7))
            return

        if self.left_p < l_x+origin_x <= self.right_p:
            right_edge = l_x+origin_x

        if l_x+origin_x > self.right_p:
            right_edge = self.right_p

        if self.left_p <= origin_x < self.right_p:
            left_edge = origin_x

        if origin_x <= self.left_p:
            left_edge =self.left_p

        N_part = int(  Nstat*((right_edge - left_edge)//dx)*(l_y//dy)  )
        lengths_p = np.array([right_edge - left_edge, l_y])
        orig_p = np.array([left_edge,origin_y])

        ## distribution of electron momentums
        v_init = np.array([v_x ,v_y])
        if EDFx=='thermal':
            velocsX = np.random.randn(N_part,1)
        elif EDFx == 'waterbag':
            velocsX = 2*(np.random.rand(N_part,1)-0.5)

        if EDFy=='thermal':
            velocsY = np.random.randn(N_part,1)
        elif EDFy == 'waterbag':
            velocsY = 2*(np.random.rand(N_part,1)-0.5)

        coords = lengths_p*np.random.rand(N_part,2) + orig_p
        velocs = v_init*np.hstack((velocsX,velocsY ))

        # defining particles IDs
        id = np.arange(N_part) + N_max*self.rank
        id.shape = (id.shape[0],1)

        gamma_part = np.ones((N_part,1))
        ewght = np.ones((N_part,1))

        self.particles = np.hstack((coords, velocs, gamma_part, ewght, id))
        self.particles = np.array(self.particles, dtype='d', order='F')
        indx = np.where(id==markID)[0]

        if tracking and indx.shape[0]!=0:		# defining initial parameters of tracked particles
            self.particles[indx,:2] = marked_coord
            self.particles[indx,2:4] = marked_veloc

    def proton_init(self):
        ''' initial distribution of protons '''
        self.particles_p = self.particles
        self.particles_p[:,2:4] = 0.
        self.particles_p[:,-1] = 1
        self.dens_ion_calcul()

    def leap_frog_motion(self):
        '''calculates the initial Lorentz-factors of the particles, relativistic density (n_e/gamma)
		   and particle velocities on the first time half-step (leap-frog scheme).'''

        # caluclating the initial electrostatic potential
        self.e_static_potent()

        # initial optical lattice potential
#        if self.particles.shape[0]!=0:
        self.pump = self.toy.pump_func(0)

        # velocity at half-step
        dummy_particles = fort_ewok.push_them_relativ( \
            e_charge*self.stat_potent, np.power(np.abs(self.pump),2), \
            self.particles[:,:5] ,self.left_p, self.bottom_p, \
            self.x_ij, self.y_ij, dx, dy, 0.5*dt)

        self.particles[:,2:5] = dummy_particles[:,2:5]
        self.particles[:,:2] = self.particles[:,:2] + dt*self.particles[:,2:4]/self.particles[:,[4]]

        part_2left, part_2right, self.particles = self.boundaries(self.particles)
        part_from_left, part_from_right = self.part_exchange(part_2left, part_2right)
        self.particles = np.vstack((self.particles, part_from_left,part_from_right))

        # calculating the electron density
        self.dens2gamma_calcul()

    def motion(self):
        ''' 1) calculates scattered and pump vector potentials and electrostatic potential to push the particles
		    2) applies sorting algorithms for left-right-bottom-top boundaries
		    3) exchanges the particles between Ewoks
		    4) calculates relativistic density (n_e/gamma)
		    1,2,4 steps are using Fortran subroutines via f2py interface generator.
        '''
        # electrostatic potential calculation
        self.e_static_potent()

        # pushing
        if self.particles.shape[0]!=0:
            self.particles[:,:5] = fort_ewok.push_them_relativ(\
                e_charge*self.stat_potent, np.abs(self.potent_em + self.pump)**2,\
                self.particles[:,:5], self.left_p, self.bottom_p,\
                self.x_ij, self.y_ij, dx, dy, dt)

        # boundary conditions
        part_2left, part_2right, self.particles = self.boundaries(self.particles)

        # exchanging particles
        part_from_left, part_from_right = self.part_exchange(part_2left, part_2right)

        # repacking particles
        self.particles = np.vstack((self.particles, part_from_left,part_from_right))

        # calculating the density divided by gamma
        self.dens2gamma_calcul()

    def motion_proton(self):
        ''' push the ions (same as electrons) '''
        if self.particles_p.shape[0]!=0:
            self.particles_p[:,:5] = fort_ewok.push_them_relativ( (-e_charge/M_p)*self.stat_potent, \
                np.abs(self.potent_em + self.pump)**2/M_p**2, \
                self.particles_p[:,:5], self.left_p, self.bottom_p, \
                self.x_ij, self.y_ij, dx, dy, dt)

        part_2left, part_2right, self.particles_p = self.boundaries(self.particles_p)
        part_from_left, part_from_right = self.part_exchange(part_2left, part_2right)
        self.particles_p = np.vstack((self.particles_p, part_from_left,part_from_right))

        self.dens_ion_calcul()

    def potential_calcul(self,time):
        ''' 1) solves the e.m. paraxial equation using explict-implicit hybrid 1-st order scheme for (x,t) and FFT for y variables
		    2) exchanges the boundary data between Ewoks
		    3) injects the seed wave
		    4) defines e.m. wave absorbtion on bottom-top boundaries for non-periodic case *y_period=False*
		    1 can be done either by NumPy or using Fortran subroutines via f2py interface generator.
        '''
        self.pump = self.toy.pump_func(time) # external ponderomotive potential (lattice)

        # transverse fourier of the right-hand-side of electromagnetic equation
        rhs_fft = fft(-1.*e_charge*self.dens2gamma*(self.pump+self.potent_em))

		# iterating electromagnetic equation
        self.potent_em_fft = fort_ewok.explicit_scheme_loop(self.potent_em_fft,self.k_y_loc,rhs_fft, dt,dt/dx)

        # right boundary potential to pass to the right
        self.potent_next[:] = +self.potent_em_fft[-1,:]
        self.comm.barrier()

        # exchange of the fields between nodes
        self.potent_exchange()

        # injecting the seed at the left boundary
        if a_seed>0.0 and self.rank == 0:
            if seed_delay<=time and time<seed_duration:
                self.potent_previous += self.toy.seed(time)

        self.potent_em_fft[0,:] = +self.potent_previous[:]

        # absorbing horizontal boundaries: field is restored to real space, multiplied by envelope, Fourier-projected
        self.potent_em = ifft(self.potent_em_fft)
        if not y_period:
            self.potent_em = self.potent_em*self.env
            self.potent_em_fft = fft(self.potent_em)

    def e_static_potent(self):
        ''' 1) gathers the density from all Ewoks to the Cheif
                    2) calculates electrostatic potential solving Poisson equation by factorization of LU operator
                    3) divides and broadcasts electrostatic potential to Ewoks
        '''
        if e_static:

            # calculate density for Poisson equation (non-relativistic)
            self.dens_calcul()

			# composing the global density distribution (mpi.gather)
            if self.rank==0:
                data_tot = np.zeros(( (self.bins[0]+1)*self.size, self.bins[1]+1), dtype='d')
            else:
                data_tot = None
            data_pass = np.array(self.data-self.dens_p, dtype='d')
            self.comm.Gather([data_pass, MPI.DOUBLE], [data_tot, MPI.DOUBLE], root=0)

			# removing the doubled columns and solving Poisson equation
            if self.rank == 0:
                data_tot = np.delete(data_tot, self.to_kick, 0)
                self.stat_potents = poiss_splu_solve(self.solver, (2*np.pi)**2*data_tot, dx, dy)
            else:
                self.stat_potents = np.zeros((self.size*self.bins[0]+1,self.bins[1]+1),dtype='d')

			# broadcasting the electrostatic potential and picking the local part
            self.comm.Bcast([self.stat_potents, MPI.DOUBLE], root=0)
            self.stat_potent = +self.stat_potents[self.rank*self.bins[0]:(self.rank+1)*self.bins[0]+1]
        else:
            self.stat_potents = np.array([], dtype='d')
            self.stat_potent = np.zeros(self.bins+1, dtype='d')

    def boundaries(self, particles):
        ''' sorts the particles for exchange returning particles going to the left, right or staying.
		    accounts for *y_particle_absorb* option for absorb particle absorbtion or periodicity at bottom-top boundaries
		    calls Fortran subroutines via f2py interface generator
        '''
        if particles.shape[0]==0:
            return np.zeros(particles.shape), np.zeros(particles.shape), np.zeros(particles.shape)
        else:
            if y_particle_absorb:
                part_2left, part_2right, part_2stay, Npart_2left, Npart_2right, Npart_2stay,  Npart_absorbed = \
                  fort_ewok.bound_migration_absorbtion(particles, self.left_p, self.right_p, self.bottom_p, self.top_p)
                self.Npart_absorbed += Npart_absorbed
            else:
                part_2left, part_2right, part_2stay, Npart_2left, Npart_2right, Npart_2stay = \
                    fort_ewok.bound_and_migration(particles, self.left_p, self.right_p, self.bottom_p, self.top_p)

            part_2left = +part_2left[:Npart_2left,:]
            part_2right = +part_2right[:Npart_2right,:]
            part_2stay = +part_2stay[:Npart_2stay,:]

            return part_2left, part_2right, part_2stay

    def wave_nums(self):
        ''' creates wave number axis for FFT methods'''
        self.k_y_full = np.fft.fftfreq(self.bins_full[1], dy)
        kx_dummy, self.k_y_loc = np.meshgrid(np.zeros(self.bins[0]+1), self.k_y_full)
        self.k_y_loc = self.k_y_loc.T

    def tansverse_env(self):
        '''creates the envelop for e.m. field absorbtion at bottom-top boundaries for *y_period=False* '''
        self.env  = np.ones(self.bins+1,dtype='d')
        pillow_x, pillow_y = np.mgrid[0:1:self.env.shape[0]*1.0j, 0:np.pi/2:pillow_size*1.0j]
        self.env[:,:pillow_size] *= np.sin(pillow_y)
        self.env[:,-pillow_size:] *= np.cos(pillow_y)
