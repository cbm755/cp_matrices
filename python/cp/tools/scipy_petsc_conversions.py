#import numpy as np
from numpy import zeros as np_zeros
#import scipy.sparse as sp
from petsc4py import PETSc
from mpi4py import MPI



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

