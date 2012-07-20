"""Example.
...
"""

import numpy as np
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import  mlab
from time import time
from operator import mod

from cp import surfaces
from cp import cpGrid, cpOps

hemisphere = surfaces.Hemisphere()

# Wrap the object in CPBar, for accurately imposing boundary
# conditions on shapes with boundaries (like a hemisphere)
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

#######

usz = M.shape[0]  # system size
Nsteps = 500      # How many time steps
dt = 0.1*dx*dx    # note O(dx^2) b/c this is explicit
kappa = 1.0       # diffusion coefficient


# setup (random) initial conditions
u0 = 1*np.random.randn(usz)
# scale to max abs 1
u0 = u0 / max(abs(u0.max()), abs(u0.min()))
u0plot = Eplot*E*u0
u = u0


# setup viz
f = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(640,640))
mlab.clf()
if True:
    # build a pipeline so we can explicitly change the data later
    src = mlab.pipeline.grid_source(x, y, z, scalars=u0plot.reshape(x.shape))
    normals = mlab.pipeline.poly_data_normals(src)
    surf = mlab.pipeline.surface(normals)
else:
    # this is easier but less convenient to change the data
    s = mlab.mesh(x, y, z, scalars=real(uplot.reshape(x.shape)))
# TODO: size is not fontsize
mlab.title('Initial conditions', size=0.2)
mlab.show()

# Explicit Forward Euler time stepping
for kt in xrange(Nsteps):
    now = time()
    utilde = E*u
    unew = utilde[0:usz] + (dt*kappa)*(D*utilde)
    print "timestep " + str(kt) + ", took %.3g s" % (time()-now)

    t = 0 + dt*(kt+1)

    # Do viz for some timesteps
    if (kt < 10) or (mod(kt+1,10) == 0):
        uplot = Eplot*(E*unew)  # brackets necessary for speed here
        # update the viz
        src.data.point_data.scalars = Eplot*(E*unew)
        src.data.point_data.scalars.name = 'scalars'
        src.data.modified()
        mlab.title("time = " + str(t) + ", step #" + str(kt+1), size=0.2)
        mlab.show()
        #mlab.savefig('hemisphere' + str(ii) + '_' + str(eval) + '.png')

        # A simpler but less efficient approach: draw a new surface each time
        if False:
            mlab.clf()
            s = mlab.mesh(x, y, z, scalars=uplot.reshape(x.shape))
            #, vmin=-absmax, vmax=absmax)
            mlab.title("time = " + str(t) + ", step #" + str(kt+1))
            mlab.show()
            pause(0.1)

    u = unew


