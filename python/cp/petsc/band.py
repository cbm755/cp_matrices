'''
Created on Jul 29, 2012

@author: nullas
'''
from __future__ import division
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
import exceptions
try:
    from cp.cpOps import buildInterpWeights
except ImportError:
    from cpOps import buildInterpWeights

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
        self.m = opt.get('m',4)
        self.StencilWidth = opt.get('sw',1)
        self.Dim = opt.get('d',3)
        self.interpDegree = opt.get('p',3)
        self.hBlock = 4/self.M
        self.hGrid = self.hBlock/self.m
        self.sizes = (self.M,)*self.Dim
        self.dx = self.hGrid

    def SelectBlock(self,surface = None):

        if surface is None:
            surface = self.surface
        comm = self.comm
        numTotalBlock = self.M**self.Dim
        numBlockAssigned = numTotalBlock // comm.size + int(comm.rank < (numTotalBlock % comm.size))
        Blockstart = comm.exscan(numBlockAssigned)
        if comm.rank == 0:Blockstart = 0
        indBlock = np.arange(Blockstart,Blockstart+numBlockAssigned)
        subBlock = self.BlockInd2SubWithoutBand(indBlock)
        BlockCenterCar = self.BlockSub2CenterCarWithoutBand(subBlock)
        cp,_,_,_ = surface.cp(BlockCenterCar)
        dBlockCenter = self.norm1(cp-BlockCenterCar)
        p = self.interpDegree
        if p % 2 == 1:
            p = ( p + 1 ) / 2
        else:
            p = ( p + 2 ) / 2
        bw = 1.1*((p+2)*self.hGrid+self.hBlock/2)#*np.sqrt(self.Dim)
        (lindBlockWithinBand,) = np.where(dBlockCenter<bw)
        lindBlockWithinBand = lindBlockWithinBand+Blockstart
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
        self.ni2pi = PETSc.AO().createMapping(self.gindBlockWBand.getArray().astype(np.int64))

    def getCoordinatesWithGhost(self):
        '''Return the coordinates of local vector.'''
        # TODO: this may return repeated superfluous coordinates.
        leng =  self.m+self.StencilWidth*2
        x = np.arange(leng**self.Dim)
        x = self.Ind2Sub(x, (leng,)*self.Dim).astype(np.double)
        x *= (self.hBlock+self.hGrid*self.StencilWidth*2)/leng
        x -= (self.StencilWidth-1/2)*self.hGrid
        x = np.tile(x,(self.numBlockWBandAssigned,1))


        y = self.gindBlockWBand.getArray()
        y = self.BlockInd2CornerCarWithoutBand(y)
        y = np.tile(y,leng**self.Dim)
        y = y.reshape((-1,self.Dim))

        return x + y

    def computeCP(self):
        cp = self.getCoordinates()
        cp,_,_,_ = self.surface.cp(cp)
        self.cp = cp
    def test_initialu(self,f):
        self.gvec.setArray(f(self.getCoordinates()))


    def initialu(self,f):
        self.gvec.setArray(f(self.cp))

    def getCoordinates(self):
        '''Return the coordinates of global vector.'''
        leng =  self.m
        x = np.arange(leng**self.Dim)
        x = self.Ind2Sub(x, (leng,)*self.Dim).astype(np.double)
        x *= (self.hBlock)/leng
        x += self.hGrid/2
        x = np.tile(x,(self.numBlockWBandAssigned,1))


        y = self.gindBlockWBand.getArray()
        y = self.BlockInd2CornerCarWithoutBand(y)
        y = np.repeat(y,leng**self.Dim,axis=0)

        return x + y

    def findIndForIntpl(self,cp):
        '''find the indices of interpolation points'''
        #find base point first
        #subBlock is the sub of Block the base points lie in
        p = self.interpDegree
        if p%2 == 1:
            offset = p // 2
        else:
            offset = p // 2 - 1
        subBlock = np.floor_divide(cp+2-self.dx/2,self.dx)
        bp = (subBlock-(offset-1/2))*self.dx - 2
        subInBlock = np.mod(subBlock,self.m)
        subBlock = np.floor_divide(subBlock,self.m)
#        corner = self.BlockSub2CornerCarWithoutBand(subBlock)

        subInBlock -= offset

        offsetBlock = np.floor_divide(subInBlock,self.m)
        subInBlock = np.mod(subInBlock,self.m)
        subBlock += offsetBlock




        p = self.interpDegree + 1
        d = self.Dim
        x = np.arange(p**d)
        x = self.Ind2Sub(x, (p,)*d)
#        x = np.tile(x,(cp.shape[0],1))
        #time consuming or memory consuming? choose one
        subInBlock = np.repeat(subInBlock,p**d,axis=0)
        subBlock = np.repeat(subBlock,p**d,axis=0)
#        offset = np.zeros(subBlock.shape)
#        x += subInBlock
        subInBlock += np.tile(x,(cp.shape[0],1))
        subBlock += np.floor_divide(subInBlock,self.m)

#        ind = np.where( x > self.m )
#        offset[ind] = 1
        subInBlock = np.mod(subInBlock,self.m)

        indBlock = self.BlockSub2IndWithoutBand(subBlock)
        indBlock = self.ni2pi.app2petsc(indBlock)
        ind = self.Sub2Ind(subInBlock, (self.m,)*d)
        ind += indBlock*self.m**d
        return bp,ind.reshape((-1,p**d))

    def createGLVectors(self):

        self.SelectBlock()
        self.computeCP()



        self.larray = np.zeros((self.numBlockWBandAssigned,)+\
                               (self.m+self.StencilWidth*2,)*self.Dim,order='F')
        self.lvec = PETSc.Vec().createWithArray(self.larray,comm=self.comm)
        lsize = self.numBlockWBandAssigned*self.m**self.Dim
        self.gvec = PETSc.Vec().createMPI((lsize,PETSc.DECIDE))
        self.gvec.setUp()
        self.wvec = self.gvec.copy()


#        self.createIndicesHelper()
        tind = np.arange((self.m+self.StencilWidth*2)**self.Dim)
        tind = tind.reshape((self.m+self.StencilWidth*2,)*self.Dim,order='F')

        for dim in xrange(self.Dim):
            tind = np.delete(tind,0,dim)
            tind = np.delete(tind,np.s_[-1],dim)

        tind = tind.flatten(order='F')

#        ISList = []
#        c = (self.m+self.StencilWidth*2)**self.Dim
#        for i in xrange(self.BlockWBandStart,self.BlockWBandEnd):
#            ti = i*c
#            ISList.extend(list(tind+ti))
        tind = np.tile(tind,self.numBlockWBandAssigned)
        ttind = np.arange(self.BlockWBandStart,self.BlockWBandEnd)
        tt = (self.m+2*self.StencilWidth)**self.Dim
        ttind *= tt
        ttind = np.repeat(ttind, self.m**self.Dim)
        ttind = tind + ttind


        ISFrom = PETSc.IS().createGeneral(ttind,comm=self.comm)
        self.l2g = PETSc.Scatter().create(self.lvec,ISFrom,self.gvec,None)

        #generate scatter global2local
        tind = np.arange(tt)
        tind = self.Ind2Sub(tind,(self.m+2*self.StencilWidth,)*self.Dim)
        tind -= self.StencilWidth
        tind = np.tile(tind,(self.numBlockWBandAssigned,1))
        ttind = self.ni2pi.petsc2app(np.arange(self.BlockWBandStart,self.BlockWBandEnd))
        ttind = self.BlockInd2SubWithoutBand(ttind)
        ttind = np.repeat(ttind,tt,axis=0)
        ttind += tind/self.m
        tind = np.mod(tind,self.m)
        tind = self.Sub2Ind(tind, (self.m,)*self.Dim)
        ttind = self.BlockSub2IndWithoutBand(ttind)
        ttind = self.ni2pi.app2petsc(ttind)
        ttind *= self.m**self.Dim
        tind += ttind
        (ind,) = np.where(tind>=0)
        tind = tind[ind]
        ISTo = ind+self.BlockWBandStart*tt
        ISFrom = PETSc.IS().createGeneral(tind)
        ISTo = PETSc.IS().createGeneral(ISTo)
        self.g2l = PETSc.Scatter().create(self.gvec,ISFrom,self.lvec,ISTo)
        return self.larray,self.lvec,self.gvec,self.wvec




    def toZero(self,gvec = None):
        if gvec is None:
            gvec = self.gvec
        tozero,zvec = PETSc.Scatter.toZero(self.gvec)# return values are self.tozero, self.zvec
        tozero.scatter(gvec, zvec, PETSc.InsertMode.INSERT)
        return zvec

    @staticmethod
    def toZeroStatic(gvec):
        tozero,zvec = PETSc.Scatter.toZero(gvec)
        tozero.scatter(gvec, zvec, PETSc.InsertMode.INSERT)
        return zvec

    def createAnyMat(self,rp,weights,NNZ = None):
        if NNZ is None:
            NNZ = (rp.shape[0],rp.shape[0])
        shape0 = rp.shape[0]
        tt = self.m**self.Dim
        start = self.BlockWBandStart
        size = (self.m,)*self.Dim
        rpt = np.tile(rp,(tt,1))

        m = PETSc.Mat().create(comm=self.comm)
        m.setSizes((self.wvec.sizes,self.gvec.sizes))
        m.setFromOptions()
        m.setPreallocationNNZ(NNZ)
        ind = np.arange(tt)
        tsubInBlock = self.Ind2Sub(ind, size)
        tsubInBlock = np.repeat(tsubInBlock,shape0,axis=0)
        tsubInBlock += rpt
        offset = np.floor_divide(tsubInBlock,self.m)
        subInBlock = np.mod(tsubInBlock,self.m)
        ones = np.ones(tt*shape0,dtype=np.int64)
        indInBlock = self.Sub2Ind(subInBlock, size)
        for block in xrange(self.numBlockWBandAssigned):
            tx = (block+start)*tt
            index = ind + tx
            #petsc ->  natural -> sub -> +offset -> natural -> petsc
            nind = self.ni2pi.petsc2app(block + start)
            nind = ones*nind
            nsub = self.BlockInd2SubWithoutBand(nind)
            nsub += offset
            nind = self.BlockSub2IndWithoutBand(nsub)
            pind = self.ni2pi.app2petsc(nind)

            pind *= tt
            pind += indInBlock
#            pind = pind.reshape((-1,shape0))
            for i in xrange(tt):
                m[index[i],pind[shape0*i:shape0*(i+1)]] = weights
        m.assemble()
        return m
    
    def createAnyMatGL(self, rp, weights, NNZ=None):
        if NNZ is None:
            NNZ = (rp.shape[0],0)
        shape0 = rp.shape[0]
        tt = self.m**self.Dim
        ltt = (self.m + 2 * self.StencilWidth) ** self.Dim
        start = self.BlockWBandStart
        size = (self.m, ) * self.Dim
        lsize = (self.m + 2 * self.StencilWidth, ) * self.Dim
        rpt = np.tile(rp,(tt,1))

        m = PETSc.Mat().create(comm = self.comm)
        m.setSizes((self.gvec.sizes,self.lvec.sizes))
        m.setFromOptions()
        m.setPreallocationNNZ(NNZ)
        
        ind = np.arange(tt)
        tsubInBlock = self.Ind2Sub(ind, size)
        tsubInBlock = np.repeat(tsubInBlock,shape0,axis=0)
        tsubInBlock += rpt
        subInBlock = tsubInBlock + self.StencilWidth
        ones = np.ones(tt*shape0,dtype=np.int64)
        indInBlock = self.Sub2Ind(subInBlock, lsize)
        
        for block in xrange(self.numBlockWBandAssigned):
            tx = (block + start) * tt
            index = ind + tx
            #petsc ->  natural -> sub -> +offset -> natural -> petsc
#            nind = self.ni2pi.petsc2app(block + start)
#            nind = ones*nind
#            nsub = self.BlockInd2SubWithoutBand(nind)
#            nind = self.BlockSub2IndWithoutBand(nsub)
            pind = ones * (block + start)

            pind *= ltt
            pind += indInBlock
#            pind = pind.reshape((-1,shape0))
            for i in xrange(tt):
                m[index[i],pind[shape0*i:shape0*(i+1)]] = weights
        m.assemble()
        return m
    
    



    def createExtensionMat(self,cp = None):
        '''create a real PETSc.Mat() for extension'''
        p = self.interpDegree + 1
        d = self.Dim

        if cp is None:
            wvec = self.wvec
            cp = self.cp
        else:
            wvec = PETSc.Vec().createMPI((cp.shape[0],PETSc.DECIDE))
        gvec = self.gvec

        extMat = PETSc.Mat().create(self.comm)
        extMat.setSizes((wvec.sizes,gvec.sizes))
        extMat.setFromOptions()
        extMat.setPreallocationNNZ((p**d,p**d))


        Xgrid,ind = self.findIndForIntpl(cp)
        weights = buildInterpWeights(Xgrid,cp,self.hGrid,p)

        (start,end) = extMat.getOwnershipRange()
        ranges = end - start


        for i in xrange(ranges):
            extMat[i+start,ind[i]] = weights[i]



        extMat.assemble()
        return extMat

    def createExtensionMatForLoop(self,cp = None):
        '''create a real PETSc.Mat() for extension using for loop'''
        p = self.interpDegree + 1
        d = self.Dim

        if cp is None:
            wvecsizes = self.wvec.sizes
            cp = self.cp
        else:
            wvecsizes = (cp.shape[0],PETSc.DECIDE)
        gvec = self.gvec

        extMat = PETSc.Mat().create(self.comm)
        extMat.setSizes((wvecsizes,gvec.sizes))
        extMat.setFromOptions()
        extMat.setPreallocationNNZ((p**d,p**d))


        (start,end) = extMat.getOwnershipRange() #@UnusedVariable

        bsize = 1000
        for i in xrange(0,cp.shape[0],bsize):
            Xgrid,ind = self.findIndForIntpl(cp[i:i+bsize])
            weights = buildInterpWeights(Xgrid,cp[i:i+bsize],self.hGrid,p)
            ranges = weights.shape[0]


            for j in xrange(ranges):
                extMat[j+start+i,ind[j]] = weights[j]



        extMat.assemble()
        return extMat



    @staticmethod
    def Ind2Sub(index,size):
        if np.isscalar(index):
            ind = np.int64(index)
        else:
            ind = index.copy()
            if ind.dtype is not np.dtype('int64'):
                ind = ind.astype(np.int64)
        total = 1
        for i in size:
            total *= i
        if np.isscalar(ind):
            if ind < 0 or ind > total:
                raise exceptions.IndexError('BlockInd2SubWithoutBand')
            return np.column_stack(np.unravel_index(ind,size,order='F'))[0]

        else:
            if ind.any() < 0 or ind.any() > total:
                raise exceptions.IndexError('BlockInd2SubWithoutBand')
            return np.column_stack(np.unravel_index(ind,size,order='F'))

    @staticmethod
    def Sub2Ind(tsub,size):
        if tsub.dtype is not np.dtype('int'):
            sub = tsub.T.astype(np.int64)
        else:
            sub = tsub.T
        return np.ravel_multi_index(sub,size,order='F')

    def BlockInd2SubWithoutBand(self, indices):
        indices = np.unravel_index(indices.astype(np.int64),
                                   (self.M,)*self.Dim, order='F')
        return np.column_stack(indices)

    def BlockSub2IndWithoutBand(self,tsub):
        if tsub.ndim == 1:
            sub = tsub.reshape((-1,self.Dim)).T
        else:
            sub = tsub.T
        if sub.dtype is not np.dtype('int'):
            sub = sub.astype(np.int64)
        return np.ravel_multi_index(sub,(self.M,)*self.Dim,order='F')


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
            return np.absolute(x).max()
        return np.amax(np.absolute(x),axis=1)






