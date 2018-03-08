import numpy as np
from mpi4py import MPI
from input import *

class Exchange:
    def potent_exchange(self):
        '''exchange of e.m. potential between Ewoks (blocking)
            accounts for *x_period* condition for e.m. field
        '''

        self.potent_previous = np.zeros(self.potent_next.size, dtype ='complex')

        if self.size==1:
            return

        if np.mod(self.rank,2)==0 and self.rank != self.size-1:
            self.comm.Send([self.potent_next.flatten(), MPI.COMPLEX], dest=self.rank+1)
        elif np.mod(self.rank,2)==1 and self.rank != 0:
            self.comm.Recv([self.potent_previous, MPI.COMPLEX], source=self.rank-1)

        if np.mod(self.rank,2)==1 and self.rank != self.size-1:
            self.comm.Send([self.potent_next.flatten(), MPI.COMPLEX], dest=self.rank+1)
        elif np.mod(self.rank,2)==0 and self.rank != 0:
            self.comm.Recv([self.potent_previous, MPI.COMPLEX], source=self.rank-1)

        if x_period == True:
            if self.rank == self.size-1:
                self.comm.Send([self.potent_next.flatten(), MPI.COMPLEX], dest=0)
            elif self.rank == 0:
                self.comm.Recv([self.potent_previous, MPI.COMPLEX], source=self.size-1)

        self.potent_previous = self.potent_previous.reshape(self.potent_next.shape)

        self.comm.barrier()

    def dens_exchange(self, dens_2right, dens_2left):
        '''exchanges data of additional cells at boundaries, used for the 3rd order interpolation (blocking)'''
        dens_from_right = np.zeros(dens_2right.flatten().shape,dtype ='d')
        dens_from_left = np.zeros(dens_2left.flatten().shape,dtype ='d')

        if self.size==1:
            return dens_from_right.reshape(dens_2left.shape), dens_from_left.reshape(dens_2right.shape)

        if np.mod(self.rank,2)==0 and self.rank != self.size-1:
            self.comm.Send([dens_2right.flatten(), MPI.DOUBLE],dest=self.rank+1)
        elif np.mod(self.rank,2)==1 and self.rank != 0:
            self.comm.Recv([dens_from_left,MPI.DOUBLE],source=self.rank-1)

        if np.mod(self.rank,2)==1 and self.rank != self.size-1:
            self.comm.Send([dens_2right.flatten(), MPI.DOUBLE],dest=self.rank+1)
        elif np.mod(self.rank,2)==0 and self.rank != 0:
            self.comm.Recv([dens_from_left,MPI.DOUBLE],source=self.rank-1)

        if np.mod(self.rank,2)==0 and self.rank != 0:
            self.comm.Send([dens_2left.flatten(), MPI.DOUBLE],dest=self.rank-1)
        elif np.mod(self.rank,2)==1 and self.rank != self.size-1:
            self.comm.Recv([dens_from_right,MPI.DOUBLE],source=self.rank+1)

        if np.mod(self.rank,2)==1 and self.rank != 0:
            self.comm.Send([dens_2left.flatten(), MPI.DOUBLE],dest=self.rank-1)
        elif np.mod(self.rank,2)==0 and self.rank != self.size-1:
            self.comm.Recv([dens_from_right,MPI.DOUBLE],source=self.rank+1)

        self.comm.barrier()
        if not x_particle_absorb:

            if self.rank == self.size-1:
                self.comm.Send([dens_2right.flatten(), MPI.DOUBLE],dest=0)
            elif self.rank == 0:
                self.comm.Recv([dens_from_left,MPI.DOUBLE],source=self.size-1)

            if self.rank == 0:
                self.comm.Send([dens_2left.flatten(), MPI.DOUBLE],dest=self.size-1)
            elif self.rank == self.size-1:
                self.comm.Recv([dens_from_right,MPI.DOUBLE],source=0)

        self.comm.barrier()

        return  dens_from_right.reshape(dens_2left.shape), dens_from_left.reshape(dens_2right.shape)

    def part_exchange(self, part_2left,part_2right):
        '''exchange of the particles between Ewoks (blocking)'''

        shape_2left = np.array(part_2left.shape)
        shape_2right = np.array(part_2right.shape)
        shape_from_left = np.array([0,part_2left.shape[1]],dtype='int')
        shape_from_right = np.array([0,part_2left.shape[1]],dtype='int')

        if self.size==1:
            return np.zeros(shape_from_left),  np.zeros(shape_from_right)

        # SHAPES EXCHAGE
        if np.mod(self.rank, 2)==0 and self.rank != self.size-1:
            self.comm.Send([shape_2right, MPI.INTEGER], dest=self.rank+1)
        elif np.mod(self.rank,2)==1 and self.rank != 0:
            self.comm.Recv([shape_from_left, MPI.INTEGER], source=self.rank-1)

        if np.mod(self.rank, 2)==1 and self.rank != self.size-1:
            self.comm.Send([shape_2right, MPI.INTEGER], dest=self.rank+1)
        elif np.mod(self.rank, 2)==0 and self.rank != 0:
            self.comm.Recv([shape_from_left, MPI.INTEGER],source=self.rank-1)

        if self.rank == self.size-1 and not x_particle_absorb:
            self.comm.Send([shape_2right,MPI.INTEGER],dest=0)
        elif self.rank == 0 and not x_particle_absorb:
            self.comm.Recv([shape_from_left,MPI.INTEGER],source=self.size-1)

        if np.mod(self.rank,2)==0 and self.rank != 0:
            self.comm.Send([shape_2left,MPI.INTEGER],dest=self.rank-1)
        elif np.mod(self.rank,2)==1 and self.rank != self.size-1:
            self.comm.Recv([shape_from_right,MPI.INTEGER],source=self.rank+1)

        if np.mod(self.rank,2)==1 and self.rank != 0:
            self.comm.Send([shape_2left,MPI.INTEGER],dest=self.rank-1)
        elif np.mod(self.rank,2)==0 and self.rank != self.size-1:
            self.comm.Recv([shape_from_right,MPI.INTEGER],source=self.rank+1)

        if self.rank == 0:
            self.comm.Send([shape_2left,MPI.INTEGER],dest=self.size-1)
        elif self.rank == self.size-1:
            self.comm.Recv([shape_from_right,MPI.INTEGER],source=0)

        if self.rank == self.size-1:
            self.comm.Send([shape_2right,MPI.INTEGER],dest=0)
        elif self.rank == 0:
            self.comm.Recv([shape_from_left,MPI.INTEGER],source=self.size-1)

        self.comm.barrier()

####################################################################################
        part_from_left = np.zeros(shape_from_left)
        part_from_right = np.zeros(shape_from_right)
###########################	PARTICLE EXCHANGE	################################

        if np.mod(self.rank,2)==0 and self.rank != self.size-1:
            self.comm.Send([part_2right,MPI.DOUBLE],dest=self.rank+1)
        elif np.mod(self.rank,2)==1 and self.rank != 0:
            self.comm.Recv([part_from_left,MPI.DOUBLE],source=self.rank-1)

        if np.mod(self.rank,2)==1 and self.rank != self.size-1:
            self.comm.Send([part_2right,MPI.DOUBLE],dest=self.rank+1)
        elif np.mod(self.rank,2)==0 and self.rank != 0:
            self.comm.Recv([part_from_left,MPI.DOUBLE],source=self.rank-1)

        if np.mod(self.rank,2)==0 and self.rank != 0 :
            self.comm.Send([part_2left,MPI.DOUBLE],dest=self.rank-1)
        elif np.mod(self.rank,2)==1 and self.rank != self.size-1:
            self.comm.Recv([part_from_right,MPI.DOUBLE],source=self.rank+1)

        if np.mod(self.rank,2)==1 and self.rank != 0 :
            self.comm.Send([part_2left,MPI.DOUBLE],dest=self.rank-1)
        elif np.mod(self.rank,2)==0 and self.rank != self.size-1:
            self.comm.Recv([part_from_right,MPI.DOUBLE],source=self.rank+1)

        self.comm.barrier()

        if x_particle_absorb:
            if self.rank == self.size-1:
                part_from_right = np.zeros((0,part_2right.shape[1]))
            elif self.rank == 0:
                part_from_left = np.zeros((0,part_2right.shape[1]))

        else:

            if self.rank == self.size-1:
                part_2right[:,0]-=(right-left)
                self.comm.Send([part_2right,MPI.DOUBLE],dest=0)
            elif self.rank == 0:
                self.comm.Recv([part_from_left,MPI.DOUBLE],source=self.size-1)

            if self.rank == 0:
                part_2left[:,0] += (right-left)
                self.comm.Send([part_2left,MPI.DOUBLE],dest=self.size-1)
            elif self.rank == self.size-1:
                self.comm.Recv([part_from_right,MPI.DOUBLE],source= 0)
            self.comm.barrier()

        return part_from_left, part_from_right




########################################################################
########################################################################
########################################################################

    def dens_exchange_orig(self, dens_2right, dens_2left):
        '''exchanges data of additional cells at boundaries, used for the 3rd order interpolation (blocking)'''
        dens_from_right = np.zeros(dens_2right.shape,dtype ='d')
        dens_from_left = np.zeros(dens_2left.shape,dtype ='d')

        if self.size==1:
            return dens_from_right, dens_from_left

        if np.mod(self.rank,2)==0 and self.rank != self.size-1:
            self.comm.Send([dens_2right, MPI.DOUBLE],dest=self.rank+1)
        elif np.mod(self.rank,2)==1 and self.rank != 0:
            self.comm.Recv([dens_from_left,MPI.DOUBLE],source=self.rank-1)

        if np.mod(self.rank,2)==1 and self.rank != self.size-1:
            self.comm.Send([dens_2right,MPI.DOUBLE],dest=self.rank+1)
        elif np.mod(self.rank,2)==0 and self.rank != 0:
            self.comm.Recv([dens_from_left,MPI.DOUBLE],source=self.rank-1)

        if np.mod(self.rank,2)==0 and self.rank != 0:
            self.comm.Send([dens_2left,MPI.DOUBLE],dest=self.rank-1)
        elif np.mod(self.rank,2)==1 and self.rank != self.size-1:
            self.comm.Recv([dens_from_right,MPI.DOUBLE],source=self.rank+1)

        if np.mod(self.rank,2)==1 and self.rank != 0:
            self.comm.Send([dens_2left,MPI.DOUBLE],dest=self.rank-1)
        elif np.mod(self.rank,2)==0 and self.rank != self.size-1:
            self.comm.Recv([dens_from_right,MPI.DOUBLE],source=self.rank+1)

        self.comm.barrier()
        if not x_particle_absorb:

            if self.rank == self.size-1:
                self.comm.Send([dens_2right,MPI.DOUBLE],dest=0)
            elif self.rank == 0:
                self.comm.Recv([dens_from_left,MPI.DOUBLE],source=self.size-1)

            if self.rank == 0:
                self.comm.Send([dens_2left,MPI.DOUBLE],dest=self.size-1)
            elif self.rank == self.size-1:
                self.comm.Recv([dens_from_right,MPI.DOUBLE],source=0)

        self.comm.barrier()

        return  dens_from_right, dens_from_left


    def potent_exchange_orig(self):
        '''exchange of e.m. potential between Ewoks (blocking)
            accounts for *x_period* condition for e.m. field
        '''

        self.potent_previous = np.zeros(np.shape(self.potent_next),dtype ='complex')

        if self.size==1:
            return

        if np.mod(self.rank,2)==0 and self.rank != self.size-1:
            self.comm.Send([self.potent_next, MPI.COMPLEX],dest=self.rank+1)
        elif np.mod(self.rank,2)==1 and self.rank != 0:
            self.comm.Recv([self.potent_previous, MPI.COMPLEX],source=self.rank-1)

        if np.mod(self.rank,2)==1 and self.rank != self.size-1:
            self.comm.Send([self.potent_next,MPI.COMPLEX],dest=self.rank+1)
        elif np.mod(self.rank,2)==0 and self.rank != 0:
            self.comm.Recv([self.potent_previous,MPI.COMPLEX],source=self.rank-1)

        if x_period == True:
            if self.rank == self.size-1:
                self.comm.Send([self.potent_next,MPI.COMPLEX],dest=0)
            elif self.rank == 0:
                self.comm.Recv([self.potent_previous,MPI.COMPLEX],source=self.size-1)

        self.comm.barrier()
