r"""
Solve the diffusion (heat) equation

     $u_t = \kappa \triangle_s u$

on a surface.  $\triangle_s$ is the Laplace--Beltrami operator.
$\kappa$ is a scalar coefficient.

Uses implicit backward Euler timestepping.
"""

from scipy.sparse import identity
from scipy.sparse.linalg import gmres
#from scipy.sparse.linalg import gmres as gmres
#from scipy.linsolve import spsolve
#import scipy.sparse.linalg.dsolve as spsolve
from scipy.sparse.linalg import spsolve
import numpy as np

try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab
from time import time
from operator import mod

usz = M.shape[0]  # system size
# Tfinal = 1
Nsteps = 100      # How many time steps
dt = 0.1*dx       # note O(dx) b/c this is implicit
kappa = 1.0       # diffusion coefficient
Iterative = True  # using  an iterative solver is faster


# setup initial conditions
u0 = 1*np.random.randn(usz)
# scale to max abs 1
u0 = u0 / max(abs(u0.max()), abs(u0.min()))
u0plot = Eplot*E*u0
u = u0


# setup viz
f = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(640,640))
mlab.clf()
# build a pipeline so we can explicitly change the data later
src = mlab.pipeline.grid_source(x, y, z, scalars=u0plot.reshape(x.shape))
normals = mlab.pipeline.poly_data_normals(src)
surf = mlab.pipeline.surface(normals)
# TODO: size is not fontsize
mlab.title('Initial conditions', size=0.2)
mlab.show()
pause(0.1)


# Playing with preconditioning
#dds = scipy.sparse.spdiags(1.0 / A_fe.diagonal(), 0, usz,usz)
#Precond = dds
# TODO: try playing with spilu as a preconditioner (need scipy 0.8)


# Build the backward Euler time-stepping matrix
A_be = identity(usz) - (kappa*dt)*M


# Implicit Backward Euler time stepping
for kt in range(0,Nsteps):
    now = time()
    if Iterative:  # Use iterative solver
        (unew,flag) = gmres(A_be, u, \
                                x0=u, tol=1e-14, restrt=10)
        #(unew,flag) = scipy.sparse.linalg.lgmres(A_be, u, \
        #                        x0=u, tol=1e-14)
        #(unew,flag) = scipy.sparse.linalg.bicgstab(A_be, u, \
        #                         x0=u, tol=1e-14)
        if (flag != 0):
            print "  flags = " + str((flag))
    else:  # Use direct solver
        unew = spsolve(A_fe, u)

    print "timestep " + str(kt) + ", took %.3g s" % (time()-now)

    t = 0 + dt*(kt+1)

    # Do viz for some timesteps
    if (kt < 10) or (mod(kt+1,20) == 0):
        uplot = Eplot*(E*unew)  # brackets necessary for speed here
        # update the viz
        src.data.point_data.scalars = Eplot*(E*unew)
        src.data.point_data.scalars.name = 'scalars'
        src.data.modified()
        mlab.title("time = " + str(t) + ", step #" + str(kt+1), size=0.2)
        mlab.show()
        #mlab.savefig('hemisphere' + str(ii) + '_' + str(eval) + '.png')

    u = unew

