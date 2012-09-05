r"""
Solve the diffusion (heat) equation

     $u_t = \kappa \triangle_s u$

on a surface.  $\triangle_s$ is the Laplace--Beltrami operator.
$\kappa$ is a scalar coefficient.

Uses implicit backward Euler timestepping.
"""

import numpy as np
from scipy.sparse import identity
from scipy.sparse.linalg import gmres, spsolve
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab
from time import time

from cp import surfaces, cpGrid, cpOps


hemisphere = surfaces.Hemisphere()

# Wrap the object in CPBar, for accurately imposing boundary
# conditions on shapes with boundaries (like a hemisphere)
# See [Macdonald, Brandman, Ruuth]
hemisphere = surfaces.CPBar(hemisphere)

# base grid point
x0 = -2. * np.ones(3)
# Grid spacing from which it'll start refining
initialdx = 4.
# Maximum levels of refinement
maxlevel = 4

TreeGrid = cpGrid.CPGrid('Spamname', hemisphere.cp, hemisphere._dim, x0,
                         initialdx, interp_degree=3, levels=maxlevel+1)


# Hemisphere parametrization in cartesian coordinates, for example for
# plotting
x, y, z = hemisphere.ParamGrid(64)
PlotPts = np.vstack((x.ravel(), y.ravel(), z.ravel())).T

TreeGrid.findStencilSets(maxlevel)
TreeGrid.buildListsFromStencilSets(maxlevel)

def boundary_function(bdy):
    if bdy == 1 or bdy == 2:
        return 'dirichlet_2nd_order'
    else:
        return None

TreeGrid.findStencilsOnStencilSets(maxlevel, boundary_function)

Lev = TreeGrid.Levolve[maxlevel]
dx = Lev[0].dx
Lex = TreeGrid.Lextend[maxlevel]
Grid = TreeGrid.Grids[maxlevel]

D = cpOps.buildDiffMatrix(Lev, Lex)
E = cpOps.buildExtensionMatrix(Lev, Lex)
Eplot = cpOps.buildEPlotMatrix(Grid, Lev, Lex, PlotPts, interp_degree=3)


# the implicit closest point method matrix.  Almost a product of D and
# E (see Macdonald and Ruuth 2009)
M = cpOps.LinearDiagonalSplitting(D, E)


usz = M.shape[0]  # system size
# Tfinal = 1
Nsteps = 100      # How many time steps
dt = 0.1*dx       # note O(dx) b/c this is implicit
kappa = 1.0       # diffusion coefficient
Iterative = True  # using  an iterative solver is much (~40x) faster


# setup (random) initial conditions
u0 = np.random.randn(usz)
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


# Playing with preconditioning
#dds = scipy.sparse.spdiags(1.0 / A_fe.diagonal(), 0, usz,usz)
#Precond = dds
# TODO: try playing with spilu as a preconditioner (need scipy 0.8)


# Build the backward Euler time-stepping matrix
A_be = identity(usz) - (kappa*dt)*M


# Implicit Backward Euler time stepping
for kt in xrange(Nsteps):
    now = time()
    if Iterative:  # Use iterative solver
        (unew,flag) = gmres(A_be, u,
                            x0=u, tol=1e-14, restrt=10)
        #(unew,flag) = scipy.sparse.linalg.lgmres(A_be, u,
        #                        x0=u, tol=1e-14)
        #(unew,flag) = scipy.sparse.linalg.bicgstab(A_be, u,
        #                         x0=u, tol=1e-14)
        if flag != 0:
            print "  flags = " + str((flag))
    else:  # Use direct solver
        unew = spsolve(A_be, u)

    print "timestep " + str(kt) + ", took %.3g s" % (time()-now)

    t = 0 + dt*(kt+1)

    # Do viz for some timesteps
    if kt < 10 or not (kt+1) % 20:
        uplot = Eplot*(E*unew)  # brackets necessary for speed here
        # update the viz
        src.data.point_data.scalars = Eplot*(E*unew)
        src.data.point_data.scalars.name = 'scalars'
        src.data.modified()
        mlab.title("time = " + str(t) + ", step #" + str(kt+1), size=0.2)
        mlab.show()
        #mlab.savefig('hemisphere' + str(ii) + '_' + str(eval) + '.png')
    u = unew

