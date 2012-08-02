'''
Created on Jul 29, 2012

@author: nullas
'''
from __future__ import division
from petsc4py import PETSc
import scipy as sp
from numpy import array as a
import exceptions


class Band(object):
    '''
    select from coarse grid.
    '''


    def __init__(self,surface,comm):
        '''
        Constructor
        '''
        self.comm = comm
        self.surface = surface
        OptDB = PETSc.Options()
        self.M = OptDB.getInt('M',20)
        self.m = OptDB.getInt('m',10)
        self.StencilWidth = OptDB.getInt('p',4)
        self.Dim = OptDB.getInt('d',3)
        self.hBlock = 4/self.M
        self.hGrid = self.hBlock/self.m
        
        numTotalBlocks = self.M**self.Dim
        numBlocksAssigned = numTotalBlocks // comm.size + int(comm.rank < (numTotalBlocks % comm.size))
        BlockStart = comm.exscan(numBlocksAssigned)
        if comm.rank == 0:BlockStart = 0
        indBlocks = sp.arange(BlockStart,BlockStart+numBlocksAssigned)
        subBlocks = self.BlockInd2SubWithoutBand(indBlocks)
        BlockCenterCar = self.BlockSub2CenterCarWithoutBand(subBlocks)
        _,dBlockCenter,_,_ = surface.cp(BlockCenterCar)
        (lindBlockWithinBand,) = sp.where(dBlockCenter<1/2*\
                                          (self.hBlock-self.hGrid)*sp.sqrt(self.Dim))
        lindBlockWithinBand = lindBlockWithinBand+BlockStart
        lBlockSize = a(lindBlockWithinBand.size,'i')
        numTotalBlocksWBand = comm.allreduce(lBlockSize)
        self.numTotalBlocksWBand = numTotalBlocksWBand
        self.numBlockWBandAssigned = numTotalBlocksWBand // comm.size + int(comm.rank < (numTotalBlocksWBand % comm.size))
        gindBlocksWBand = PETSc.Vec().createMPI((lBlockSize,PETSc.DECIDE))
        gindBlocksWBand.setArray(lindBlockWithinBand)
        self.gindBlockWBand = PETSc.Vec().createMPI((self.numBlockWBandAssigned,PETSc.DECIDE))
        ISAll = PETSc.IS().createStride(numTotalBlocksWBand,step=1,comm=comm)
        scatter = PETSc.Scatter().create(gindBlocksWBand,ISAll,self.gindBlockWBand,None) 
        scatter.scatter(gindBlocksWBand,self.gindBlockWBand,PETSc.InsertMode.INSERT)
#        PETSc.Sys.syncPrint('Process {0} got {1} Blocks'.format(comm.rank,lBlockSize))
#        PETSc.Sys.syncFlush()
#        PETSc.Sys.Print('Total Blocks {0}'.format(self.numTotalBlocksWBand))
        
    
    def BlockInd2SubWithoutBand(self,index):
        if sp.isscalar(index):
            ind = index
        else:
            ind = index.copy()
        if sp.isscalar(ind):
            if ind < 0 or ind > self.M**self.Dim-1:
                raise exceptions.IndexError('BlockInd2SubWithoutBand')
            sub = sp.zeros(self.Dim)
            for i in range(0,self.Dim):
                sub[i] = ind%self.M
                ind //= self.M
            return sub
            
        else:
            if ind.any() < 0 or ind.any() >self.M**self.Dim-1:
                raise exceptions.IndexError('BlockInd2SubWithoutBand')
            sub = sp.zeros((ind.size,self.Dim))
            for i in range(0,self.Dim):
                sub[:,i] = ind%self.M
                ind //= self.M
            return sub
    def BlockSub2CenterCarWithoutBand(self,sub):       
        return sub/self.M*4-2+self.hBlock/2
        
        
        