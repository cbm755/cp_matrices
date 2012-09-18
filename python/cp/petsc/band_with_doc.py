'''
Created on Jul 29, 2012

@author: nullas (modified by yujia)
'''
from __future__ import division
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
import exceptions
try:
    from cp.cpOps import findGridInterpBasePt, buildInterpWeights
except ImportError:
    from cpOps import findGridInterpBasePt, buildInterpWeights

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
        self.ll = opt.get('ll',(-2,)*self.Dim)

    def SelectBlock(self,surface = None):
        '''
	Starting from structured blocks in a big cube, find and keep those blocks which are near the surface.
	Before selection, need to distribute blocks to processors s.t. each processor has roughly the same
	number of blocks; and after selecting those blocks within band, some processors might have few (or no) 
	blocks but others might have many, we let PETSc decide how to re-distribute the blocks to ensure
	load balance again.
	Technically the re-distribution is done via the 'scatter' function.
	'''
        if surface is None:
            surface = self.surface
        comm = self.comm
        numTotalBlock = self.M**self.Dim
        
	# Decide how many blocks should be assigned to each processor.
	# If numTotalBlock = k*comm.size + r: each of the first r processors has k+1 blocks,
        # each of the rest processors has k blocks.
        numBlockAssigned = numTotalBlock // comm.size + int(comm.rank < (numTotalBlock % comm.size))
        
	# Exclusive scan with default operation '+'. 
        # On the i-th processor, sum up the 'numBlockAssinged' in the previous 0,1...,(i-1)-th processors.
        Blockstart = comm.exscan(numBlockAssigned)
        if comm.rank == 0:Blockstart = 0
		
	# Linear indices of the Blocks(consecutive numbers), each processor owns its own part of global indices.
        indBlock = np.arange(Blockstart,Blockstart+numBlockAssigned)
        
	# Sub indices (i,j,k) of the Blocks with respect to the virtual shape M*M*M
        subBlock = self.BlockInd2SubWithoutBand(indBlock)
        
	# Find the coordinates in the center of each block 
        BlockCenterCar = self.BlockSub2CenterCarWithoutBand(subBlock)

        cp,_,_,_ = surface.closest_point(BlockCenterCar)
        #cp,_,_,_ = surface.cp(BlockCenterCar)
        dBlockCenter = self.norm1(cp-BlockCenterCar)
        p = self.interpDegree
        if p % 2 == 1:
            p = ( p + 1 ) / 2
        else:
            p = ( p + 2 ) / 2
        bw = 1.1*((p+2)*self.hGrid+self.hBlock/2)#*np.sqrt(self.Dim)
        (lindBlockWithinBand,) = np.where(dBlockCenter<bw)

        # The above np.where finds the indices starting from 0 (local indices), we add Blockstart to make it global.  
        lindBlockWithinBand = lindBlockWithinBand+Blockstart

        lBlockSize = lindBlockWithinBand.size
        numTotalBlockWBand = comm.allreduce(lBlockSize)

        numBlockWBandAssigned = numTotalBlockWBand // comm.size + int(comm.rank < (numTotalBlockWBand % comm.size))

        # Creat a PETSc vector from numpy array? Not sure what's exactly doing here. 
	# 'lindBlockWBandFrom' is the vector FROM which we want to scatter. 
        lindBlockWBandFrom = PETSc.Vec().createWithArray(lindBlockWithinBand,comm=comm)

        # Allocate memory for the vector 'self.gindBlockWBand' TO which we want to scatter. 
        self.gindBlockWBand = PETSc.Vec().createMPI((numBlockWBandAssigned,PETSc.DECIDE),comm=comm)


#        gsubBlockWBandFrom = PETSc.Vec().createMPI((self.Dim*lBlockize,PETSc.DECIDE),comm=comm)
#        gsubBlockWBandFrom.setArray(lsubBlockWBand)
#        self.gsubBlockWBand = PETSc.Vec().createMPI((self.Dim*self.numBlockWBandAssigned,PETSc.DECIDE),comm=comm)

        BlockWBandStart = comm.exscan(numBlockWBandAssigned)
        if comm.rank == 0:
            BlockWBandStart = 0
        self.BlockWBandStart = BlockWBandStart
        
        # Index sets of the vector FROM which we want to scatter.  
        LInd = PETSc.IS().createStride(numBlockWBandAssigned,\
                                       first=BlockWBandStart,\
                                       step=1,comm=comm)
        
        # Re-scatter the blocks within band. The last argument of PETSc.Scatter().create() should be the index sets 
        # of the vector TO which we want to scatter. 'None' means that we fill the entire vector 'self.gindBlockWBand'. 
        scatter = PETSc.Scatter().create(lindBlockWBandFrom,LInd,self.gindBlockWBand,None)
        scatter.scatter(lindBlockWBandFrom,self.gindBlockWBand,PETSc.InsertMode.INSERT)
       
	# Natural order Index To Petsc order Index
        # TODO: colin needs to change this int64 to int32: what is correct thing to do here?
        self.ni2pi = PETSc.AO().createMapping(self.gindBlockWBand.getArray().astype(np.int64))
        #self.ni2pi = PETSc.AO().createMapping(self.gindBlockWBand.getArray().astype(np.int32))
		
        self.numTotalBlockWBand = numTotalBlockWBand
	self.numBlockWBandAssigned = numBlockWBandAssigned
	self.BlockWBandEnd = BlockWBandStart + numBlockWBandAssigned
		
	
    def getCoordinates(self):
	''' 
	Return the coordinates of the fine grid. 
        Each grid centered at each element(square in 2D, cube in 3D) of the Cartesian mesh.
	'''
        # find the coordinates of the fine grid in a single reference block with lower-left corner at (0,0,0).
	leng = self.m
        x = np.arange(leng**self.Dim)
        x = self.Ind2Sub(x, (leng,)*self.Dim).astype(np.double)
        x *= self.hGrid
        x += self.hGrid/2 
        x = np.tile(x,(self.numBlockWBandAssigned,1))

        # find the coordinates of all the find grid by shifting the reference coordinates.
        y = self.gindBlockWBand.getArray()
        y = self.BlockInd2CornerCarWithoutBand(y)
        y = np.repeat(y,leng**self.Dim,axis=0)

        return x + y
    
    def computeCP(self):
	cp = self.getCoordinates()
        cp,_,_,_ = self.surface.closest_point(cp)
        #cp,_,_,_ = self.surface.cp(cp)
        self.cp = cp

    def createGLVectors(self):
        '''
	Create global and local vectors for the values at fine grid.
	Currently only global vectors 'self.gvec' and 'self.wvec' is used.
	'''

        self.SelectBlock()
        
        # Create global vectors 'self.gvec' and 'self.wvec'.
        lsize = self.numBlockWBandAssigned*self.m**self.Dim
        self.gvec = PETSc.Vec().createMPI((lsize,PETSc.DECIDE))
        self.gvec.setUp()
        self.wvec = self.gvec.copy()

	# The rest lines of code for this def, except the last return line, is creating
        # local vectors, which is not used currently. Could skip reading them for now.

        self.larray = np.zeros((self.numBlockWBandAssigned,)+\
                               (self.m+self.StencilWidth*2,)*self.Dim,order='F')
        self.lvec = PETSc.Vec().createWithArray(self.larray,comm=self.comm)

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

    def test_initialu(self,f):
        self.gvec.setArray(f(self.getCoordinates()))

    def initialu(self,f):
        self.gvec.setArray(f(self.cp))
    
    def subToPetscInd(self,sub):
        '''Convert sub-indices to Petsc indices. See comments in 'createExtMat'. We need this helper function
        because we only know the PETSc indices of big blocks. So we first find a point is in which block, then
        find the natural linear index inside that block.'''
        m = self.m
        M = self.M
        dim = self.Dim
        
        # sub-indices of the Block containing the point.
        subBlock = np.floor_divide(sub,np.array((m,)*dim)[...,np.newaxis,np.newaxis])
        # sub-indices of the point in the block containing it.
        subInBlock = np.mod(sub,np.array((m,)*dim)[...,np.newaxis,np.newaxis])
        
        indBlock = np.ravel_multi_index(subBlock,(M,)*dim,order='F')
        petscIndBlock = self.ni2pi.app2petsc(indBlock)
        indInBlock = np.ravel_multi_index(subInBlock,(m,)*dim,order='F')
        ind = indInBlock + petscIndBlock*m**dim        

        return ind
        
    def createMat(self, size, ind, weights, NNZ = None):
        ''' 
        Create a PETSc matrix, ind is the global PETSc indices. 
        size: size of the matrix 
        ind: Suppose each row of the matrix has nnz non-zero elements, and this processor has k rows, 
             then 'ind' should be a k*nnz array indicating the global column indices of matrix entries.
        weights: 'weights' stores the values of matrix entries corresponding to 'ind'.
        NNZ: an array with 2 entries. NNZ[0] is the maximal number of nonzeros in each row of the whole matrix;
             NNZ[1] is the maximal number of nonzeros which might be in another processor in one row. 
        ''' 
        if NNZ is None:
            NNZ = (ind.shape[1],ind.shape[1])

        m = PETSc.Mat().create(comm=self.comm)
        m.setSizes(size)
        m.setFromOptions()
        m.setPreallocationNNZ(NNZ)
        (start,end) = m.getOwnershipRange()
        ranges = end - start

        # We specify part of the matrix rows on this processor, but the row and column indices should be global.
        for i in xrange(ranges):
            m[i+start,ind[i]] = weights[i]

        m.assemble()
        return m
    
    def createExtensionMat(self, p = None, cp = None):
        '''create a PETSc.Mat() for the closest point extension'''

        if p is None: p = self.interpDegree

        if cp is None:
            wvecsizes = self.wvec.sizes
            cp = self.cp
        else:
            wvecsizes = (cp.shape[0],PETSc.DECIDE)
        gvec = self.gvec

        dim = self.Dim
        dx = (self.hGrid,)*dim
        M = self.M
        m = self.m
        # The lower left grid point 'll' is at the center of the lower left element 'self.ll'
        ll = np.array(self.ll) + self.hGrid/2;

        # find the sub-indices of the base points, 'subBasept' is a N*dim array (N:number of base-points).
        subBasept = findGridInterpBasePt(cp, dx, ll, p)
        # offsets is the sub-indices of the interpolation stencil, a dim*STENCIL array (STENCIL=(p+1)^dim).
        # If in 'C' order, should be the following line: 
        offsets = np.mgrid[(slice(p+1),)*dim].reshape((dim, -1))
        # or equivalently the following two lines:
#        x = np.arange((p+1)**dim)
#        offsets = np.column_stack(np.unravel_index(x,(p+1,)*dim,order='C')).T
        #If in the structure branch, it should be 'Fortran' order:
#        x = np.arange((p+1)**dim)
#        offsets = np.column_stack(np.unravel_index(x,(p+1,)*dim,order='F')).T
        # The sub indices of the whole interpolation stencil with Fortran order, a dim*N*STENCIL array.         
        sub = ( subBasept[...,np.newaxis] + offsets[np.newaxis,...] ).transpose((1,0,2))
        # Convert sub indices to global PETSc indices
        ind = self.subToPetscInd(sub)

        basept = subBasept*dx + ll
        weights = buildInterpWeights(basept,cp,dx,p+1)

        E = self.createMat((wvecsizes,gvec.sizes),ind,weights)
        
        return E

    def createLaplacianMat(self):
        '''create a PETSc.Mat() for discrete Laplacian.
        TODO: more than dim 3; more than second order; different dx,dy,dz'''
        dim = self.Dim

        dx = (self.hGrid,)*dim
        ll = np.array(self.ll) + self.hGrid/2
        pt = self.getCoordinates()
        # find the sub index of itself 
        subItself = findGridInterpBasePt(pt,dx,ll,0)

        if dim == 2:
            offsets = np.array([[0,0,1,0,-1],
                                [0,1,0,-1,0]])
            weight = np.array([-4.,1.,1.,1.,1.]) / self.hGrid**2
        elif dim == 3:
            offsets = np.array([[0,0,0,1,0,0,-1],
                                [0,0,1,0,0,-1,0],
                                [0,1,0,0,-1,0,0]]) 
            weight = np.array([-6.,1.,1.,1.,1.,1.,1.]) / self.hGrid**2

        sub = ( subItself[...,np.newaxis] + offsets[np.newaxis,...] ).transpose((1,0,2))
        # Convert sub indices to global PETSc indices
        ind = self.subToPetscInd(sub)

        weights = np.tile(weight,(ind.shape[0],1)) 

        L = self.createMat((self.gvec.sizes,self.gvec.sizes),ind,weights,(2*dim+1,dim))
        
        return L

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


    '''
    Following are some defs which might be useful in the future.
    '''
    
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

    def findIndForIntpl(self,cp):
        '''find the indices of interpolation points'''
        dim = self.Dim
        dx = (self.hGrid,)*dim
        p = self.interpDegree
        M = self.M
        m = self.m
        ll = np.array(self.ll) + self.hGrid/2;
        #find base point first
        subBasept = findGridInterpBasePt(cp, dx, ll, p)
        offsets = np.mgrid[(slice(p+1),)*dim].reshape((dim, -1))
#        x = np.arange((p+1)**dim)
#        offsets = np.column_stack(np.unravel_index(x,(p+1,)*dim,order='F')).T
        sub = ( subBasept[...,np.newaxis] + offsets[np.newaxis,...] ).transpose((1,0,2))
        ind = self.subToPetscInd(sub)

        basept = subBasept*dx + ll
        return basept,ind

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


    def createAnyMatGL(self, rp, weights, NNZ=None):
	'''
	Create matrix for double-band. 
	G,L is the notation for the two bands in the implicit CPM paper Macdonald & Ruuth
	'''
	if NNZ is None:NNZ = (rp.shape[0],0)
			
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
#           nind = self.ni2pi.petsc2app(block + start)
#           nind = ones*nind
#           nsub = self.BlockInd2SubWithoutBand(nind)
#           nind = self.BlockSub2IndWithoutBand(nsub)
	    pind = ones * (block + start)
		
	    pind *= ltt
	    pind += indInBlock
#           pind = pind.reshape((-1,shape0))
	for i in xrange(tt):
            m[index[i],pind[shape0*i:shape0*(i+1)]] = weights
				
	m.assemble()
	return m
    
        
