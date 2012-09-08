#import numpy as np
from numpy import zeros as np_zeros
#import scipy.sparse as sp
from petsc4py import PETSc
from mpi4py import MPI


def save_scipy_to_petsc_ondisk(A, fname):

    #import sys, petsc4py
    #petsc4py.init(sys.argv)
    from petsc4py import PETSc

    from scipy.sparse import find

    # TODO: OO way to break A into its list component form?
    Ai, Aj, As = find(A)
    # not quite same?
    #Ai = A.indptr
    #Aj = A.indices
    #As = A.data
    nnz = len(Aj)

    #B = PETSc.Mat().create(PETSc.COMM_WORLD)
    B = PETSc.Mat().createAIJ((A.shape))
    #B.setSizes([m*n, m*n])
    #B.setSizes((A.shape))
    B.setFromOptions()
    Istart, Iend = B.getOwnershipRange()
    print (Istart, Iend)

    #A.setValue(0, 0, 3) # Insert a single value into matrix.
    #A.setValues([0, 1], [2, 3], [1, 1, 1, 1]) # Insert a 2x2 block of values into the matrix.

    # No!, but I'm sure there is a form like this
    #B.setValues(Ai, Aj, As)

    for i in xrange(0, nnz):
        I = Ai[i]
        J = Aj[i]
        val = As[i]
        #print (i,I,J,val)
        # TODO: this is silly, but I'll probably run in serial anyway...
        if (I >= Istart) and (I < Iend):
        #B.setValue(I,J,val)
            B[I,J] = val
            #print (i,I,J,val)
        else:
            print ("fail", i,I,J,val)

    B.assemblyBegin()
    B.assemblyEnd()

    x, y = B.getVecs()
    x.set(1)
    B.mult(x,y)

    viewer = PETSc.Viewer().createBinary(fname, 'w')
    viewer(B)
    #viewer = PETSc.Viewer().createBinary('vector-x.dat', 'w')
    #viewer(x)
    #viewer = PETSc.Viewer().createBinary('vector-y.dat', 'w')
    #viewer(y)

    # load
    viewer = PETSc.Viewer().createBinary(fname, 'r')
    B2 = PETSc.Mat().load(viewer)
    #viewer = PETSc.Viewer().createBinary('vector-x.dat', 'r')
    # = PETSc.Vec().load(viewer)
    #viewer = PETSc.Viewer().createBinary('vector-y.dat', 'r')
    #v = PETSc.Vec().load(viewer)

    # check
    assert B2.equal(B)
    #assert x.equal(u)
    #assert y.equal(v)



# Most of the rest of this is based on code I found at:
# http://code.google.com/p/femstab/source/browse/src/python/parfemstab.py

def csrmatrix2PETScMat(L):
    """
    Converts a sequential scipy sparse matrix (on process 0) to a PETSc
    Mat ('aij') matrix distributed on all processes
    input : L, scipy sparse matrix on proc 0
    output: PETSc matrix distributed on all procs
    """

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Get the data from the sequential scipy matrix
    if rank == 0:
        if L.format == 'csr':
            L2 = L
        else:
            L2 = L.tocsr()
        Ai  = L2.indptr
        Aj  = L2.indices
        Av  = L2.data
        nnz = len(Aj)
        n,m = L2.shape
    else:
        n   = None
        m   = None
        nnz = None
        Ai  = None
        Aj  = None
        Av  = None

    # Broadcast sizes
    n   = comm.bcast(n  ,root = 0)
    m   = comm.bcast(m  ,root = 0)
    nnz = comm.bcast(nnz,root = 0)

    B = PETSc.Mat()
    B.create(comm)
    B.setSizes([n, m])
    B.setType('aij')
    B.setFromOptions()

    # Create a vector to get the local sizes, so that preallocation can be done later
    V = PETSc.Vec()
    V.create(comm)
    V.setSizes(n)
    V.setFromOptions()
    istart,iend = V.getOwnershipRange()
    V.destroy()

    nloc = iend - istart

    Istart = comm.gather(istart,root = 0)
    Iend   = comm.gather(iend  ,root = 0)

    if rank == 0:
        nnzloc = np_zeros(comm.size,'int')
        for i in range(comm.size):
            j0        = Ai[Istart[i]]
            j1        = Ai[Iend  [i]]
            nnzloc[i] = j1 - j0
    else:
        nnzloc = None

    nnzloc = comm.scatter(nnzloc,root = 0)

    ai = np_zeros(nloc+1   ,PETSc.IntType)
    aj = np_zeros(nnzloc+1 ,PETSc.IntType)
    av = np_zeros(nnzloc+1 ,PETSc.ScalarType)

    if rank == 0:
        j0        = Ai[Istart[0]]
        j1        = Ai[Iend  [0]]
        ai[:nloc  ] = Ai[:nloc]
        aj[:nnzloc] = Aj[j0:j1]
        av[:nnzloc] = Av[j0:j1]

    for iproc in range(1,comm.size):
        if rank == 0:
            i0        = Istart[iproc]
            i1        = Iend  [iproc]
            j0        = Ai[i0]
            j1        = Ai[i1]
            comm.Send(Ai[i0:i1], dest=iproc, tag=77)
            comm.Send(Aj[j0:j1], dest=iproc, tag=78)
            comm.Send(Av[j0:j1], dest=iproc, tag=79)
        elif rank == iproc:
            comm.Recv(ai[:nloc  ], source=0, tag=77)
            comm.Recv(aj[:nnzloc], source=0, tag=78)
            comm.Recv(av[:nnzloc], source=0, tag=79)

    ai = ai- ai[0]
    ai[-1] = nnzloc+1

    B.setPreallocationCSR((ai,aj))
    B.setValuesCSR(ai,aj,av)
    B.assemble()

    return B


def array2PETScVec(v):
    """
    Converts (copies) a sequential array/vector on process 0
    to a distributed PETSc Vec
    input : v, numpy array on proc 0, None (or whatever) on other proc
    output: PETSc Vec distributed on all procs
    """

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # v is (probably) only redefined on proc 0
    if rank == 0:
        n = len(v)
    else:
        n = None

    n = comm.bcast(n, root = 0)
    #print "DEBUG", __name__, "rank=", rank, "n=", n

    x = PETSc.Vec()
    x.create(comm)
    x.setSizes(n)
    x.setFromOptions()
    istart,iend = x.getOwnershipRange()

    nloc = iend - istart
    Istart = comm.gather(istart,root = 0)
    Iend   = comm.gather(iend  ,root = 0)

    vloc = np_zeros(nloc,PETSc.ScalarType)

    if rank == 0:
        vloc[:nloc  ] = v[:nloc]

    for iproc in range(1,comm.size):
        if rank == 0:
            i0        = Istart[iproc]
            i1        = Iend  [iproc]
            comm.Send(v[i0:i1], dest=iproc, tag=77)
        elif rank == iproc:
            comm.Recv(vloc, source=0, tag=77)

    x.setArray(vloc)

    return x


def PETScVec2array(x):
    """
    Converts (copies) a distributed PETSc Vec to a sequential array on process 0
    input : x, PETSc Vec distributed on all procs
    output: numpy array on proc 0
    """

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    vloc = x.getArray()
    n    = x.getSize()

    istart,iend = x.getOwnershipRange()

    nloc = iend - istart
    Istart = comm.gather(istart,root = 0)
    Iend   = comm.gather(iend  ,root = 0)

    if rank == 0:
        v = np_zeros(n,PETSc.ScalarType)
    else:
        v = None

    if rank == 0:
        v[:nloc  ] = vloc

    for iproc in range(1,comm.size):
        if rank == 0:
            i0        = Istart[iproc]
            i1        = Iend  [iproc]
            comm.Recv(v[i0:i1], source=iproc, tag=77)
        elif rank == iproc:
            comm.Send(vloc, dest=0, tag=77)

    return v

