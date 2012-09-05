"""Example.
...
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#from matplotlib.mlab import griddata

from time import time

from cp import surfaces
from cp import cpGrid, cpOps

circle = surfaces.Circle()

# Wrap the object in CPBar, for accurately imposing boundary
# conditions on shapes with boundaries (like a hemisphere)
# See [Macdonald, Brandman, Ruuth]
# Not for a circle!
#circle = surfaces.CPBar(circle)

# base grid point
x0 = -2. * np.ones(2)
# Grid spacing from which it'll start refining
initialdx = 4.
# Maximum levels of refinement
maxlevel = 4

TreeGrid = cpGrid.CPGrid('Spamname', circle.closest_point, circle.dim, x0,
                         initialdx, interp_degree=3, levels=maxlevel+1)


# Circle parametrization in cartesian coordinates, for example for
# plotting
x, y = circle.parametric_grid(64)
PlotPts = np.vstack((x.ravel(), y.ravel(),)).T

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
# We don't use this here
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
fig = plt.figure()
ax = fig.add_subplot(111,  aspect='equal')

# s is the marker size
scatter = ax.scatter(x, y, c=u0plot, cmap=cm.jet, s=80)
ax.set_title('Initial conditions')
fig.show()
raw_input('Press sth to continue')

# Explicit Forward Euler time stepping
for kt in xrange(Nsteps):
    now = time()
    utilde = E*u
    unew = utilde[0:usz] + (dt*kappa)*(D*utilde)
    print "timestep {0: 6d} took  {1:.3f} s".format(kt, time()-now)

    t = 0 + dt*(kt+1)

    # Do viz for some timesteps
    if kt < 10 or not (kt+1) % 5:
        uplot = Eplot * (E*unew)  # brackets necessary for speed here
        scatter.set_array(uplot)
        ax.set_title("time = {0: 6.3f}, step #{1: 5d}".format(t, kt+1))
        fig.canvas.draw()
    u = unew
