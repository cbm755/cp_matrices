# Colin's petsc test.

import numpy as np
import pickle
import timeit
import sys
from mpi4py import MPI
import petsc4py
from petsc4py import PETSc
import cp.tools.scipy_petsc_conversions as conv


petsc4py.init(sys.argv)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# this will load the picke on each processor
#(dx,initial_u,final_u) = pickle.load(file('non_petsc_data.pickle'))

# alternatively, load it once and share the scalar around
if rank == 0:
    (dx, cpm, Tf, numtimesteps, dt, initial_u, final_u) = pickle.load(file('non_petsc_data.pickle'))
else:
    (dx, cpm, Tf, numtimesteps, dt, initial_u, final_u) = (None, None, None, None, None, None, None)

# Broadcast variables to other processors
(dx, cpm, Tf, numtimesteps, dt) = comm.bcast((dx, cpm, Tf, numtimesteps, dt), root = 0)

#buf = np.zeros(2, PETSc.ScalarType)
#for iproc in xrange(1,comm.size):
#    if rank == 0:
#        comm.Send(dx, dest=iproc, tag=123)
#    elif rank == iproc:
#        comm.Recv(buf, source=0, tag=123)


# which matrices should we load
if cpm == 0:
    viewer = PETSc.Viewer().createBinary('Lmatrix.dat', 'r')
    L_Mat = PETSc.Mat().load(viewer)
    viewer = PETSc.Viewer().createBinary('Ematrix.dat', 'r')
    E_Mat = PETSc.Mat().load(viewer)
elif cpm == 1:
    viewer = PETSc.Viewer().createBinary('Mmatrix.dat', 'r')
    M_Mat = PETSc.Mat().load(viewer)
elif cpm == 2:
    viewer = PETSc.Viewer().createBinary('Amatrix.dat', 'r')
    A_Mat = PETSc.Mat().load(viewer)


# this only sets the local part
#v.setArray(initial_u.copy())
# Have to do this carefully:
v = conv.array2PETScVec(initial_u)
v2 = v.copy()
#v2 = L_Mat.getVecRight()


# pre-multiply the matrices by dt.  I don't like this much.
if cpm == 0:
    # replace laplacian with dt*L
    L_Mat.scale(dt)
elif cpm == 1:
    # replace laplacian with dt*M
    M_Mat.scale(dt)


start_time = timeit.default_timer()

# various CPM algorithms
if cpm == 0:  # explicit Euler, Ruuth--Merriman
    for kt in xrange(numtimesteps):
        L_Mat.multAdd(v, v, v2)    # v2 = v + (dt*L)*v
        E_Mat.mult(v2, v)          # v = E*v2
        t = kt * dt

elif cpm == 1:  # explicit Euler, von Glehn--Maerz--Macdonald
    for kt in xrange(numtimesteps):
        M_Mat.multAdd(v, v, v2)    # v2 = v + (dt*M)*v
        # todo: how to just assign data in v2 to v?
        v2.swap(v)
        t = kt * dt

elif cpm == 2:  # implicit Euler, vGMM
    for kt in xrange(numtimesteps):
        #v2 = v + dt * M * v2
        #A*v2 = v
        #v = v2
        raise NotImplementedError('learn ksp/ts')
        t = kt*dt;


print "Times, rank=", rank, "time=", timeit.default_timer() - start_time


final_u2 = conv.PETScVec2array(v)
if rank == 0:
    maxdiff = max(abs(final_u - final_u2))
    print 'max diff in serial/parallel: ', maxdiff
    if maxdiff < 1e-13:
        print '\n**PASSED**: serial/parallel codes differ by < 1e-13\n'
    else:
        print '\n**FAILED**: serial/parallel code output differs\n'

