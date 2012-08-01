"""Example.
...
"""

import numpy as np
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import  mlab
from time import time

from cp import surfaces
from cp import cpGrid, cpOps

sphere = surfaces.Sphere()

# base grid point
x0 = -2. * np.ones(3)
# Grid spacing from which it'll start refining
initialdx = 4.
# Maximum levels of refinement
maxlevel = 4

TreeGrid = cpGrid.CPGrid('Spamname', sphere.cp, sphere._dim, x0,
                         initialdx, interp_degree=3, levels=maxlevel+1)


# Sphere parametrization in cartesian coordinates, for example for
# plotting
x, y, z = sphere.ParamGrid(64)
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
#import sys
#sys.exit()
# the implicit closest point method matrix.  Almost a product of D and
# E (see Macdonald and Ruuth 2009)
# We don't need this here
# M = cpOps.LinearDiagonalSplitting(D, E)

#######

#usz = M.shape[0]  # system size
usz = D.shape[0]
Nsteps = 500      # How many time steps
dt = 0.1*dx*dx    # note O(dx^2) b/c this is explicit
kappa = 1.0       # diffusion coefficient


# setup (random) initial conditions
u0 = np.random.randn(usz)
# scale to max abs 1
u0 = u0 / max(abs(u0.max()), abs(u0.min()))
u0plot = Eplot * (E * u0)  # Parentheses important for speed
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
raw_input("Press enter to continue...")

# Explicit Forward Euler time stepping
for kt in xrange(Nsteps):
    now = time()
    utilde = E*u
    unew = utilde[0:usz] + (dt*kappa)*(D*utilde)
    print "timestep " + str(kt) + ", took %.3g s" % (time()-now)

    t = 0 + dt*(kt+1)

    # Do viz for some timesteps
    if kt < 10 or not (kt+1) % 10:
        uplot = Eplot * (E*unew)  # brackets necessary for speed here
        # update the viz
        src.data.point_data.scalars = uplot
        src.data.point_data.scalars.name = 'scalars'
        src.data.modified()
        mlab.title("time = " + str(t) + ", step #" + str(kt+1), size=0.2)
        mlab.show()
        #mlab.savefig('sphere' + str(ii) + '_' + str(eval) + '.png')
    u = unew
