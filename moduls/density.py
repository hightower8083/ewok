import numpy as np
import fort_ewok
from input import *
from exchange import Exchange

class density_rout(Exchange):

    def dens_calcul(self):
        if self.particles.shape[0]!=0:
            data_ext = fort_ewok.density_2x(self.particles[:,:2], self.particles[:,-2], \
                self.left_p, self.bottom_p, dx, dy, self.bins[0], self.bins[1])
            self.data = np.array(+data_ext[2:-2,2:-2],dtype='d' )
            dens_2right, dens_2left =  np.array(data_ext[-3:,2:-2],dtype='d') , np.array(data_ext[:3,2:-2],dtype='d')
        else:
            self.data = np.zeros(self.bins+1, dtype='d')
            dens_2right, dens_2left = np.zeros((3,self.bins[1]+1),dtype='d'), np.zeros((3,self.bins[1]+1),dtype='d')

        dens_from_right, dens_from_left = self.dens_exchange(dens_2right, dens_2left)

        self.data[ :3,:] += dens_from_left
        self.data[-3:,:] += dens_from_right
        self.comm.barrier()

    def dens_ion_calcul(self):
        if self.particles_p.shape[0]!=0:
            data_ext = fort_ewok.density_2x(self.particles_p[:,:2], self.particles_p[:,-2], \
                self.left_p, self.bottom_p, dx, dy, self.bins[0], self.bins[1])
            self.dens_p = np.array(+data_ext[2:-2,2:-2],dtype='d' )
            dens_2right, dens_2left =  np.array(data_ext[-3:,2:-2],dtype='d'), np.array(data_ext[:3,2:-2],dtype='d')
        else:
            self.dens_p = np.zeros(self.bins+1,dtype='d')
            dens_2right, dens_2left = np.zeros((3,self.bins[1]+1),dtype='d'), np.zeros((3,self.bins[1]+1),dtype='d')

        dens_from_right, dens_from_left = self.dens_exchange(dens_2right, dens_2left)

        self.dens_p[:3,:] +=  dens_from_left
        self.dens_p[-3:,:] +=  dens_from_right
#       self.dens_p = self.max_wght*self.dens_p
        self.comm.barrier()

    def dens2gamma_calcul(self):
        ''' 1) calculates the 'relativistic' particles density (thanx to Rachel for the scheme)
            2) calls routine to exchange additional cells at boundaries
            uses Fortran subroutines via f2py interface generator
        '''
        if self.particles.shape[0]!=0:
            data_ext = fort_ewok.density_2x(self.particles[:,:2], \
                self.particles[:,-2]/self.particles[:,4], \
                self.left_p, self.bottom_p, dx, dy, self.bins[0], self.bins[1])

            self.dens2gamma = np.array(+data_ext[2:-2,2:-2],dtype='d' )
            dens_2right, dens_2left =  np.array(data_ext[-3:,2:-2],dtype='d') , np.array(data_ext[:3,2:-2],dtype='d')
        else:
            self.dens2gamma = np.zeros(self.bins+1, dtype='d')
            dens_2right, dens_2left = np.zeros((3,self.bins[1]+1),dtype='d'), np.zeros((3,self.bins[1]+1),dtype='d')

        dens_from_right, dens_from_left = self.dens_exchange(dens_2right, dens_2left)

        self.dens2gamma[ :3,:] += dens_from_left
        self.dens2gamma[-3:,:] += dens_from_right
        self.comm.barrier()


##########################################################################################
##########################################################################################
##########################################################################################


    def dens2gamma_calcul_tmp(self):
        ''' 1) calculates the 'relativistic' particles density (thanx to Rachel for the scheme)
            2) calls routine to exchange additional cells at boundaries
            uses Fortran subroutines via f2py interface generator
        '''
        if self.particles.shape[0]!=0:
            data_ext = fort_ewok.density_2x(self.particles[:,:2], self.particles[:,-2]/self.particles[:,4], \
                self.left_p, self.bottom_p,dx,dy, self.bins[0], self.bins[1])
            self.dens2gamma = np.array( data_ext[2:-2,2:-2], dtype='d' )
            dens_2right, dens_2left =  np.array(data_ext[-3:,2:-2], dtype='d') , np.array(data_ext[:3,2:-2],dtype='d')
        else:
            self.dens2gamma = np.zeros(self.bins+1,dtype='d')
            dens_2right, dens_2left = np.zeros((3,self.bins[1]+1), dtype='d'), np.zeros((3,self.bins[1]+1), dtype='d')

        dens_from_right, dens_from_left = self.dens_exchange(dens_2right, dens_2left)

        self.dens2gamma[:3,:] +=  dens_from_left
        self.dens2gamma[-3:,:] +=  dens_from_right
#       self.dens2gamma = self.max_wght*self.dens2gamma
        self.comm.barrier()

    def dens_calcul_tmp(self):
        ''' 1) calculates the particles density (thanx to Rachel for the sub-routine)
                    2) calls routine to exchange additional cells at boundaries
                    uses Fortran subroutines via f2py interface generator
        '''

        if self.particles.shape[0]!=0:
            data_ext = fort_ewok.density_2x(self.particles[:,:2], self.particles[:,-2], \
                self.left_p, self.bottom_p, dx, dy, self.bins[0], self.bins[1])
            self.data = np.array( data_ext[2:-2,2:-2], dtype='d' )
            dens_2right, dens_2left =  np.array(data_ext[-3:,2:-2],dtype='d'), np.array(data_ext[:3,2:-2],dtype='d')
        else:
            self.data = np.zeros(self.bins+1,dtype='d')
            dens_2right, dens_2left = np.zeros((3,self.bins[1]+1),dtype='d'), np.zeros((3,self.bins[1]+1),dtype='d')

        dens_from_right, dens_from_left = self.dens_exchange(dens_2right, dens_2left)

        self.data[:3,:] +=  dens_from_left
        self.data[-3:,:] +=  dens_from_right
#       self.data = self.max_wght*self.data
        self.comm.barrier()

