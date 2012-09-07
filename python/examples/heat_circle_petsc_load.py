# Colin's petsc test.  Currently only works on one processor: n=1
# because of mistakes in distributing the soln vector

import numpy as np
#from matplotlib import pylab as pl
from mpi4py import MPI
#from cp.surfaces.MeshWrapper import MeshWrapper
#from cp.surfaces.Sphere import Sphere
#from cp.petsc.band import Band
import sys
import petsc4py
from petsc4py import PETSc
#import matplotlib.tri as tri
import pickle

petsc4py.init(sys.argv)

(dx,initial_u,final_u) = pickle.load(file('non_petsc_data.pickle'))


viewer = PETSc.Viewer().createBinary('Lmatrix.dat', 'r')
L_Mat = PETSc.Mat().load(viewer)

viewer = PETSc.Viewer().createBinary('Ematrix.dat', 'r')
E_Mat = PETSc.Mat().load(viewer)


v = L_Mat.getVecRight()
v2 = L_Mat.getVecRight()
#v3 = L_Mat.getVecRight()




Tf = 2
dt = 0.1 * np.min(dx)**2
numtimesteps = int(Tf // dt + 1)

# TODO: only sets the local part
v.setArray(initial_u.copy())
# replace laplacian with dt*L
L_Mat.scale(dt)

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

#final_err_petsc = max_error

# TODO: only gets the local part
maxdiff = max(abs(final_u - v.getArray()))
print('max diff in serial/parallel: ', maxdiff)


