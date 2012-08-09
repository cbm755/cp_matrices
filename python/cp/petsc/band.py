'''
Created on Jul 29, 2012

@author: nullas
'''
from __future__ import division
from mpi4py import MPI
from petsc4py import PETSc
import scipy as sp
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


    def __init__(self,surface,comm,opt = None):
        '''
        Constructor
        '''
        if opt is None:
            opt = {}
        self.comm = comm
        self.surface = surface
        self.M = opt.get('M',20)
        self.m = opt.get('m',2)
        self.StencilWidth = opt.get('sw',1)
        self.Dim = opt.get('d',3)
        self.interpDegree = opt.get('p',4)
        self.hBlock = 4/self.M
        self.hGrid = self.hBlock/self.m
        self.sizes = (self.M,)*self.Dim
        
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
        p = self.interpDegree
        if p % 2 == 1:
            p = ( p + 1 ) / 2
        else:
            p = ( p + 2 ) / 2
        bw = 1.0001*(p*self.hGrid+(self.hBlock-self.hGrid)/2)*sp.sqrt(self.Dim)
        (lindBlockWithinBand,) = sp.where(dBlockCenter<bw)
        lindBlockWithinBand = lindBlockWithinBand+Blocktart
        lBlockSize = lindBlockWithinBand.size
        numTotalBlockWBand = comm.allreduce(lBlockSize)
        
        

        numBlockWBandAssigned = numTotalBlockWBand // comm.size + int(comm.rank < (numTotalBlockWBand % comm.size))
        
        lindBlockWBandFrom = PETSc.Vec().createWithArray(lindBlockWithinBand,comm=comm)
        self.gindBlockWBand = PETSc.Vec().createMPI((numBlockWBandAssigned,PETSc.DECIDE),comm=comm)
        
        

        

        
        
#        gsubBlockWBandFrom = PETSc.Vec().createMPI((self.Dim*lBlockize,PETSc.DECIDE),comm=comm)
#        gsubBlockWBandFrom.setArray(lsubBlockWBand)
#        self.gsubBlockWBand = PETSc.Vec().createMPI((self.Dim*self.numBlockWBandAssigned,PETSc.DECIDE),comm=comm)
        BlockWBandStart = comm.exscan(numBlockWBandAssigned)
        if comm.rank == 0:
            BlockWBandStart = 0
        self.BlockWBandStart = BlockWBandStart
        
                
        LInd = PETSc.IS().createStride(numBlockWBandAssigned,\
                                       first=BlockWBandStart,\
                                       step=1,comm=comm)


        self.numTotalBlockWBand = numTotalBlockWBand
        self.numBlockWBandAssigned = numBlockWBandAssigned
        
        BlockWBandEnd = BlockWBandStart + numBlockWBandAssigned
        self.BlockWBandEnd = BlockWBandEnd
        
        
        scatter = PETSc.Scatter().create(lindBlockWBandFrom,LInd,self.gindBlockWBand,None) 
        scatter.scatter(lindBlockWBandFrom,self.gindBlockWBand,PETSc.InsertMode.INSERT)
        #Natural order Index To Petsc order Index
        self.ni2pi = PETSc.AO().createMapping(self.gindBlockWBand.getArray().astype(sp.int64))
        
    def getCoordinatesWithGhost(self):
        '''Return the coordinates of local vector.'''
        # TODO: this may return repeated superfluous coordinates.
        leng =  self.m+self.StencilWidth*2
        x = sp.arange(leng**self.Dim)
        x = self.Ind2Sub(x, (leng,)*self.Dim).astype(sp.double)
        x *= (self.hBlock+self.hGrid*self.StencilWidth*2)/leng
        x -= (self.StencilWidth-1/2)*self.hGrid
        x = sp.tile(x,(self.numBlockWBandAssigned,1))
        
        
        y = self.gindBlockWBand.getArray()
        y = self.BlockInd2CornerCarWithoutBand(y)
        y = sp.tile(y,leng**self.Dim)
        y = y.reshape((-1,self.Dim))
        
        return x + y
    def getCoordinates(self):
        '''Return the coordinates of global vector.'''
        leng =  self.m
        x = sp.arange(leng**self.Dim)
        x = self.Ind2Sub(x, (leng,)*self.Dim).astype(sp.double)
        x *= (self.hBlock)/leng
        x -= self.hGrid/2
        x = sp.tile(x,(self.numBlockWBandAssigned,1))
        
        
        y = self.gindBlockWBand.getArray()
        y = self.BlockInd2CornerCarWithoutBand(y)
        y = sp.repeat(y,leng**self.Dim,axis=0)
        
        return x + y
        
    def findIndForIntpl(self,cp):
        '''find the indices of interpolation points'''
        #find base point first
        #subBlock is the sub of Block the base points lie in
        subBlock = sp.floor_divide(cp+2,self.hBlock)
        corner = self.BlockSub2CornerCarWithoutBand(subBlock)
        p = self.interpDegree
        if p%2 == 1:
            offset = p // 2
        else:
            offset = p // 2 - 1
        subInBlock = sp.floor_divide(cp-corner-self.hGrid/2,self.hGrid)
        subInBlock -= offset
        offsetBlock = sp.floor_divide(subInBlock,self.m)
        subInBlock = sp.mod(subInBlock,self.m)
        subBlock += offsetBlock
        
        
        p = self.interpDegree + 1
        d = self.Dim
        x = sp.arange(p**d)
        x = self.Ind2Sub(x, (p,)*d)
        x = sp.tile(x,(cp.shape[0],1))
        #time consuming or memory consuming? choose one
        subInBlock = sp.tile(subInBlock,p**d)
        subInBlock = subInBlock.reshape((-1,d))
        subBlock = sp.tile(subBlock,p**d)
        subBlock = subBlock.rehshape((-1,d))
        offset = sp.zeros(subBlock.shape)
        x += subInBlock
        ind = sp.where( x > self.m )
        offset[ind] = 1
        subInBlock = sp.mod(x,self.m)
        subBlock += offset
        indBlock = self.BlockSub2IndWithoutBand(subBlock)
        try:
            indBlock = self.ni2pi.app2petsc(indBlock)
        except Exception:
            print 'Not every block used for interpolation is selected\n'
        ind = self.Sub2Ind(subInBlock, (self.m,)*d)
        ind += indBlock*self.m**d
        return ind.shape((-1,p**d))        
    @staticmethod    
    def interp1d(bp,cp,N):
        from scipy.misc import comb
        w = [comb(N,i)*(-1)**i for i in range(N+1)]
        w = sp.array(w)
        
        
        
        
        
        
              
    def createGLVectors(self):
        
        self.SelectBlock()
        

        
        self.larray = sp.zeros((self.numBlockWBandAssigned,)+\
                               (self.m+self.StencilWidth*2,)*self.Dim,order='F')
        self.lvec = PETSc.Vec().createWithArray(self.larray,comm=self.comm)
        lsize = self.numBlockWBandAssigned*self.m**self.Dim
        self.gvec = PETSc.Vec().createMPI((lsize,PETSc.DECIDE))
        self.gvec.setUp()
        self.wvec = self.gvec.copy()
        
        
#        self.createIndicesHelper()
        tind = sp.arange((self.m+self.StencilWidth*2)**self.Dim)
        tind = tind.reshape((self.m+self.StencilWidth*2,)*self.Dim,order='F')
        
        for dim in xrange(self.Dim):
            tind = sp.delete(tind,0,dim)
            tind = sp.delete(tind,sp.s_[-1],dim)
            
        tind = tind.flatten(order='F')
        
#        ISList = []
#        c = (self.m+self.StencilWidth*2)**self.Dim
#        for i in xrange(self.BlockWBandStart,self.BlockWBandEnd):
#            ti = i*c
#            ISList.extend(list(tind+ti))
        tind = sp.tile(tind,self.numBlockWBandAssigned)
        ttind = sp.arange(self.BlockWBandStart,self.BlockWBandEnd)
        tt = (self.m+2*self.StencilWidth)**self.Dim
        ttind *= tt
        ttind = sp.repeat(ttind, self.m**self.Dim)
        ttind = tind + ttind
        
            
        ISFrom = PETSc.IS().createGeneral(ttind,comm=self.comm)
        self.l2g = PETSc.Scatter().create(self.lvec,ISFrom,self.gvec,None)
        
        #generate scatter global2local
        tind = sp.arange(tt)
        tind = self.Ind2Sub(tind,(self.m+2*self.StencilWidth,)*self.Dim)
        tind -= self.StencilWidth
        tind = sp.tile(tind,(self.numBlockWBandAssigned,1))
        ttind = self.ni2pi.petsc2app(sp.arange(self.BlockWBandStart,self.BlockWBandEnd))
        ttind = self.BlockInd2SubWithoutBand(ttind)
        ttind = sp.repeat(ttind,tt,axis=0)
        ttind += tind/self.m
        tind = sp.mod(tind,self.m)
        tind = self.Sub2Ind(tind, (self.m,)*self.Dim)
        ttind = self.BlockSub2IndWithoutBand(ttind)
        ttind = self.ni2pi.app2petsc(ttind)
        ttind *= tt
        tind += ttind
        ind,_ = sp.where(tind>=0)
        ISTo = ind+self.BlockWBandStart*tt
        ISFrom = PETSc.IS().createGeneral(tind)
        ISTo = PETSc.IS().createGeneral(ISTo)
        self.g2l = PETSc.Scatter(self.gvec,ISFrom,self.lvec,ISTo)
        
        
        
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
        

    @staticmethod
    def Ind2Sub(index,size):
        if sp.isscalar(index):
            ind = sp.int64(index)
        else:
            ind = index.copy()
            if ind.dtype is not sp.dtype('int64'):
                ind = ind.astype(sp.int64)
        total = 1
        for i in size:
            total *= i
        if sp.isscalar(ind):
            if ind < 0 or ind > total:
                raise exceptions.IndexError('BlockInd2SubWithoutBand')
            return sp.column_stack(sp.unravel_index(ind,size,order='F'))[0]
            
        else:
            if ind.any() < 0 or ind.any() > total:
                raise exceptions.IndexError('BlockInd2SubWithoutBand')
            return sp.column_stack(sp.unravel_index(ind,size,order='F'))
        
    @staticmethod
    def Sub2Ind(sub,size): 
        return sp.column_stack(sp.ravel_multi_index(sub.T,size,order='F'))     
                        
                
    
    def BlockInd2SubWithoutBand(self,index):
        if sp.isscalar(index):
            ind = index
        else:
            ind = index.copy()
            if ind.dtype is not sp.dtype('int64'):
                ind = ind.astype(sp.int64)
        if sp.isscalar(ind):
            if ind < 0 or ind > self.M**self.Dim-1:
                raise exceptions.IndexError('BlockInd2SubWithoutBand')
            return sp.column_stack(sp.unravel_index(ind,(self.M,)*self.Dim,order='F'))[0]
            
        else:
            if ind.any() < 0 or ind.any() > self.M**self.Dim-1:
                raise exceptions.IndexError('BlockInd2SubWithoutBand')
            return sp.column_stack(sp.unravel_index(ind,(self.M,)*self.Dim,order='F'))
        
    def BlockSub2IndWithoutBand(self,tsub):
        if tsub.ndim == 1:
            sub = tsub.reshape((-1,self.Dim)).T
        else:
            sub = tsub.T
        return sp.ravel_multi_index(sub,(self.M,)*self.Dim,order='F')
        
        
    def BlockSub2CenterCarWithoutBand(self,sub):       
        return sub/self.M*4-2+self.hBlock/2
    
    def BlockSub2CornerCarWithoutBand(self,sub):
        return sub/self.M*4-2
    
    def BlockInd2CenterCarWithoutBand(self,ind):
        return self.BlockSub2CenterCarWithoutBand(self.BlockInd2SubWithoutBand(ind))
    
    def BlockInd2CornerCarWithoutBand(self,ind):
        return self.BlockSub2CornerCarWithoutBand(self.BlockInd2SubWithoutBand(ind))
    
    
    @staticmethod
    def norm1(x):
        if x.ndim == 1:
            return sp.absolute(x).max()
        return sp.amax(sp.absolute(x),axis=1)
    
    
        
        
        
        