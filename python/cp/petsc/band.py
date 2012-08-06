'''
Created on Jul 29, 2012

@author: nullas
'''
from __future__ import division
from mpi4py import MPI
from petsc4py import PETSc
import scipy as sp
from numpy import array 
import exceptions


#class harray(array):
#    '''Providing hash for Array. Use this for indexing dictionary.'''
#    def __init__(self,obj, dtype=None, copy=True, order=None, subok=False, ndmin=0):
#        self.a = a.__init__(self,obj, dtype, copy, order, subok, ndmin)
#        
#    def __hash__(self):
#        return hashlib.sha1(self.a)
        

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
        self.m = OptDB.getInt('m',2)
        self.StencilWidth = OptDB.getInt('p',1)
        self.Dim = OptDB.getInt('d',3)
        self.hBlock = 4/self.M
        self.hGrid = self.hBlock/self.m
        
    def SelectBlock(self,surface = None):
        
        if surface is None:
            surface = self.surface
        comm = self.comm
        numTotalBlock = self.M**self.Dim
        numBlockAssigned = numTotalBlock // comm.size + int(comm.rank < (numTotalBlock % comm.size))
        Blocktart = comm.exscan(numBlockAssigned)
        if comm.rank == 0:Blocktart = 0
        indBlock = sp.arange(Blocktart,Blocktart+numBlockAssigned)
        subBlock = self.BlockInd2SubWithoutBand(indBlock)
        BlockCenterCar = self.BlockSub2CenterCarWithoutBand(subBlock)
        _,dBlockCenter,_,_ = surface.cp(BlockCenterCar)
        (lindBlockWithinBand,) = sp.where(dBlockCenter<1/2*\
                                          (self.hBlock-self.hGrid)*sp.sqrt(self.Dim))
        lindBlockWithinBand = lindBlockWithinBand+Blocktart
        lBlockSize = lindBlockWithinBand.size
        lsubBlockWBand = self.BlockInd2SubWithoutBand(lindBlockWithinBand)
        numTotalBlockWBand = comm.allreduce(lBlockSize)
        
        

        numBlockWBandAssigned = numTotalBlockWBand // comm.size + int(comm.rank < (numTotalBlockWBand % comm.size))
               
        lsubBlockWBandFrom = PETSc.Vec().createWithArray(lsubBlockWBand,comm=comm)
        self.gsubBlockWBand = PETSc.Vec().createMPI((self.Dim*numBlockWBandAssigned,PETSc.DECIDE),comm=comm)
        
        

        

        
        
#        gsubBlockWBandFrom = PETSc.Vec().createMPI((self.Dim*lBlockize,PETSc.DECIDE),comm=comm)
#        gsubBlockWBandFrom.setArray(lsubBlockWBand)
#        self.gsubBlockWBand = PETSc.Vec().createMPI((self.Dim*self.numBlockWBandAssigned,PETSc.DECIDE),comm=comm)
        BlockWBandStart = comm.exscan(numBlockWBandAssigned)
        if comm.rank == 0:
            BlockWBandStart = 0
        self.BlockWBandStart = BlockWBandStart
        
                
        LInd = PETSc.IS().createStride(numBlockWBandAssigned*self.Dim,\
                                       first=BlockWBandStart*self.Dim,\
                                       step=1,comm=comm)


        self.numTotalBlockWBand = numTotalBlockWBand
        self.numBlockWBandAssigned = numBlockWBandAssigned
        
        BlockWBandEnd = BlockWBandStart + numBlockWBandAssigned
        self.BlockWBandEnd = BlockWBandEnd
        scatter = PETSc.Scatter().create(lsubBlockWBandFrom,LInd,self.gsubBlockWBand,None) 
        scatter.scatter(lsubBlockWBandFrom,self.gsubBlockWBand,PETSc.InsertMode.INSERT)
        
        
    def createGLVectors(self):
        
        self.SelectBlock()
        

        
        self.larray = sp.zeros((self.numBlockWBandAssigned,)+\
                               (self.m+self.StencilWidth*2,)*self.Dim,order='F')
        self.lvec = PETSc.Vec().createWithArray(self.larray,comm=self.comm)
        lsize = self.numBlockWBandAssigned*self.m**self.Dim
        self.gvec = PETSc.Vec().createMPI((lsize,PETSc.DECIDE))
        self.gvec.setUp()
        
#        self.createIndicesHelper()
        tind = sp.arange((self.m+self.StencilWidth*2)**self.Dim)
        tind = tind.reshape((self.m+self.StencilWidth*2,)*self.Dim,order='F')
        
        for dim in xrange(self.Dim):
            tind = sp.delete(tind,0,dim)
            tind = sp.delete(tind,sp.s_[-1],dim)
            
        tind = tind.flatten(order='F')
        
        ISList = []
        c = (self.m+self.StencilWidth*2)**self.Dim
        for i in xrange(self.BlockWBandStart,self.BlockWBandEnd):
            ti = i*c
            ISList.extend(list(tind+ti))
            
        ISFrom = PETSc.IS().createGeneral(ISList,comm=self.comm)
        self.l2g = PETSc.Scatter().create(self.lvec,ISFrom,self.gvec,None)
            
        
        
        
        
        
#        self.gsubBlockWBand.setArray(self.asubBlockWBan[self.Dim*BlockWBandStart:self.Dim*BlockWBandEnd])
        
        
#        gindBlockWBandFrom = PETSc.Vec().createMPI((lBlockize,PETSc.DECIDE))
#        gindBlockWBandFrom.setArray(lindBlockWithinBand)
#        self.gindBlockWBand = PETSc.Vec().createMPI((self.numBlockWBandAssigned,PETSc.DECIDE))

#        PETSc.Sys.syncPrint('Process {0} got {1} Block'.format(comm.rank,lBlockize))
#        PETSc.Sys.syncFlush()
#        PETSc.Sys.Print('Total Block {0}'.format(self.numTotalBlockWBand))


#    def createIndicesHelper(self):
#        toall,vsubBlockWBand = PETSc.Scatter().toAll(self.gsubBlockWBand)
#        toall.scatter(self.gsubBlockWBand,vsubBlockWBand)
#        asubBlockWBand = vsubBlockWBand.getArray().reshape(-1,self.Dim).astype(sp.int64)
#        numTotalBlockWBand = self.numTotalBlockWBand
#        self.sub2indWBand = {tuple(asubBlockWBand[i]):i for i in xrange(numTotalBlockWBand)}
#        self.ind2subWBand = {i:asubBlockWBand[i] for i in xrange(numTotalBlockWBand)}
#        x = sp.arange(3**self.Dim)
#        x = sp.vstack(sp.unravel_index(x,(3,)*self.Dim,order='F')).T
#        x -= sp.ones(self.Dim,dtype=sp.int64)
#        lBlockSize = self.numBlockWBandAssigned
#        for i in xrange(lBlockSize):
#            sub = self.ind2subWBand[i]
#            for offset in x:
#                tsub = tuple(sub + offset)
                
    def createExtensionMat(self):     
        '''create a real PETSc.Mat() for extension'''       
        extMat = PETSc.Mat().create()
        extMat.setO
        extMat.setSizes((self.lvec.sizes,self.gvec.sizes))
        
        
                        
                
    
    def BlockInd2SubWithoutBand(self,index):
        if sp.isscalar(index):
            ind = index
        else:
            ind = index.copy()
        if sp.isscalar(ind):
            if ind < 0 or ind > self.M**self.Dim-1:
                raise exceptions.IndexError('BlockInd2SubWithoutBand')
            return sp.column_stack(sp.unravel_index(ind,(self.M,)*self.Dim,order='F'))[0]
            
        else:
            if ind.any() < 0 or ind.any() > self.M**self.Dim-1:
                raise exceptions.IndexError('BlockInd2SubWithoutBand')
            return sp.column_stack(sp.unravel_index(ind,(self.M,)*self.Dim,order='F'))
        
    def BlockSub2CenterCarWithoutBand(self,sub):       
        return sub/self.M*4-2+self.hBlock/2
        
        
        