# Colin's petsc test.

import numpy as np
import pickle
import timeit
import sys
from mpi4py import MPI
import petsc4py
from petsc4py import PETSc
import cp.tools.scipy_petsc_conversions as conv


#from cp.surfaces.MeshWrapper import MeshWrapper
#from cp.surfaces.Sphere import Sphere
#from cp.petsc.band import Band
#import matplotlib.tri as tri


petsc4py.init(sys.argv)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


# this will load the picke on each processor
#(dx,initial_u,final_u) = pickle.load(file('non_petsc_data.pickle'))

# alternatively, load it once and share the scalar around
if rank == 0:
    (dx,initial_u,final_u) = pickle.load(file('non_petsc_data.pickle'))
else:
    initial_u = None
    final_u = None
    dx = None

# Broadcast dx to other processors
dx = comm.bcast(dx, root = 0)

#buf = np.zeros(2, PETSc.ScalarType)
#for iproc in xrange(1,comm.size):
#    if rank == 0:
#        comm.Send(dx, dest=iproc, tag=123)
#    elif rank == iproc:
#        comm.Recv(buf, source=0, tag=123)

viewer = PETSc.Viewer().createBinary('Lmatrix.dat', 'r')
L_Mat = PETSc.Mat().load(viewer)

viewer = PETSc.Viewer().createBinary('Ematrix.dat', 'r')
E_Mat = PETSc.Mat().load(viewer)

# NO! only sets the local part
#v.setArray(initial_u.copy())
# Have to do this carefully:
v = conv.array2PETScVec(initial_u)

#v = L_Mat.getVecRight()
v2 = L_Mat.getVecRight()
#v3 = L_Mat.getVecRight()



Tf = 2
dt = 0.1 * np.min(dx)**2
numtimesteps = int(Tf // dt + 1)




# replace laplacian with dt*L
L_Mat.scale(dt)

start_time = timeit.default_timer()

for kt in xrange(numtimesteps):
    # v2 = v + (dt*L)*v
    L_Mat.multAdd(v, v, v2)

    E_Mat.mult(v2, v)  # v = E*v2

    t = kt * dt
    #if not kt%100 or kt == (numtimesteps-1):
    #    uplot = Eplot * v.getArray()
    #    true_solution = np.exp(-t) * np.cos(th_plot + np.pi / 2)
    #    max_error = (np.abs(true_solution - uplot)).max()
    #    print "time: {0:2f} ({1:2.2f}%), err={2:g}".format(t, 100*float(kt) / numtimesteps, max_error)


print "Times, rank=", rank, "time=", timeit.default_timer() - start_time


# TODO: only gets the local part
final_u2 = conv.PETScVec2array(v)
if rank == 0:
    maxdiff = max(abs(final_u - final_u2))
    print('max diff in serial/parallel: ', maxdiff)


