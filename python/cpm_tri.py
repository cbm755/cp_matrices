import numpy

from numpy import array as a
import cpGrid
reload(cpGrid)
import cpOps
reload(cpOps)
import stencils

import pylab
from pylab import plot

#from scipy.linalg import norm
#import scipy.sparse
#import scipy.linsolve
#import scipy.linalg

from time import time

#global SURF_SR
# SURF_SR = 1.0
# import cp_sphere
# cp_sphere.init(center=a([0,0,0]), radius=SURF_SR)
# cpfun_inner = cp_sphere.cp_sphere

# import cp_hemisphere
# cp_hemisphere.init(center=a([0,0,0]), radius=SURF_SR)
# cpfun_inner = cp_hemisphere.cp_hemisphere

import closestpoint

#q2 = closestpoint.Circle()

#q2 = closestpoint.Semicircle()
#q2 = closestpoint.ParamCurve()
#q2 = closestpoint.CircularArc()
#q2 = closestpoint.SplineCurve()
q2 = closestpoint.TriangulationSlow()

# wrap the object in CPBar, for accurately imposing boundary
# conditions on shapes with boundaries (like a hemisphere)
#q = closestpoint.CPBar(parent=q2)
q = q2

cpfun = q.cp


# This function that helps with boundary conditions when the surface
# has a boundary.  TODO: ties into "CPBar" boundary conditions as in
# Macdonald, Brandman, Ruuth preprint.  TODO: think about cleaning
# this up.
def bdyfcn(bdy):
    if bdy==1 or bdy==2:
        #return 'dirichlet_1st_order'
        #return 'dirichlet_2nd_order'
        return 'neumann_1st_order'
        #return 'neumann_2nd_order'
    else:
        return None

print 'start'


refx = numpy.array([-4.0, -4.0, -4.0])
dx = 0.5

f96 = numpy.float96
#x = numpy.array([f96(-4.0), f96(-4.0)])
#initialdx = f96(16.0)

#x = numpy.array([-2.0, -2.0])
#initialdx = 8.0

dim = len(x)
#maxlev = 8

diffstencil = stencils.Laplacian_2nd; interpdeg = 3
#diffstencil = stencils.Laplacian_4th; interpdeg = 5

Grid = cpGrid.CPFlatGrid('bunny', cpfun, dim, x, dx, interp_degree=interpdeg, diffInfo=diffstencil)



#TreeGrid = cpGrid.CPGrid('test', cpfun, dim, x, initialdx, interp_degree=interpdeg, levels=2, diffInfo=diffstencil)


## a parameterized grid for plotting
#x,y,z = q.ParamGrid(rez=64)
x,y = q.ParamGrid(rez=512)

# this can be used to view like this:
#from enthought.mayavi import mlab
#s = mlab.mesh(x, y, z, scalars=z**2)
#mlab.show()

PlotPts = numpy.vstack((x,y)).transpose()


j = maxlev - 0
TreeGrid.findStencilSets(j)
TreeGrid.buildListsFromStencilSets(j)
TreeGrid.findStencilsOnStencilSets(j, bdyfcn)
dx = TreeGrid.Levolve[j][0].dx
Lev = TreeGrid.Levolve[j]
Lex = TreeGrid.Lextend[j]
Grid = TreeGrid.Grids[j]
# TODO: I think it would be better to keep the dx**2 out of this matrix
D = cpOps.buildDiffMatrix(Lev, Lex)
D2 = cpOps.buildDiffMatrixNoDX(Lev, Lex)
#D3xb,D3xf,D3yb,D3yf = g.buildDiffMatrixTempDxDyTest(j)
E = cpOps.buildExtensionMatrix(Lev, Lex)
Eplot = cpOps.buildEPlotMatrix(Grid, Lev, Lex, PlotPts, interp_degree=interpdeg)


# the implicit closest point method matrix.  Almost a product of D and
# E (see Macdonald and Ruuth 2009)
now = time()
M = cpOps.LinearDiagonalSplitting(D, E)
print "splitting took %.3g s" % (time()-now)

Mt = (D*E).astype(numpy.float64)

if (type(dx) == numpy.float96):
    # keep a copy of the original extended precision matrix
    M_f96 = M
    M = M_f96.astype(numpy.float64)



# minnz = 1e42
# M = E3
# for i in range(0,M.shape[0]):
#     for j in range(0,M.shape[1]):
#         if M[i,j] != 0:
#             minnz = min(minnz, abs(M[i,j]))
# print minnz



def pause(howlong=None):
    import time as timemod
    if howlong == None:
        raw_input('Paused, press enter to continue')
    else:
        timemod.sleep(howlong)
