"""
Makes a grid surrounding a closest point representation and builds the
implicit Closest Point Method matrix operators on it.
"""

import numpy as np

from cp import cpGrid
reload(cpGrid)
from cp import cpOps
reload(cpOps)

#import pylab
#from pylab import plot

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

from cp import surfaces

q2 = surfaces.Hemisphere()
# wrap the object in CPBar, for accurately imposing boundary
# conditions on shapes with boundaries (like a hemisphere)
q = surfaces.CPBar(q2)
cpfun = q.cp

#cpfun = cpbar
#cpfun = cpfun_inner
#cpfun = cp_sphere.cp_sphere

# This function that helps with boundary conditions when the surface
# has a boundary.  TODO: ties into "CPBar" boundary conditions as in
# Macdonald, Brandman, Ruuth preprint.  TODO: think about cleaning
# this up.
def bdyfcn(bdy):
    if bdy==1 or bdy==2:
        #return 'dirichlet_1st_order'
        return 'dirichlet_2nd_order'
        #return 'neumann_1st_order'
        #return 'neumann_2nd_order'
    else:
        return None

print 'start'


#x = numpy.array([-2.0, -2.0])
#dx = 4.0

#f96 = numpy.float96
#x = numpy.array([f96(-2.0), f96(-2.0), f96(-2.0)])
#dx = f96(8.0)

x = np.array([-2.0, -2.0, -2.0])
initialdx = 4.0

dim = len(x)
maxlev = 4
TreeGrid = cpGrid.CPGrid('test', cpfun, dim, x, initialdx, interp_degree=3, levels=maxlev+1)


## a parameterized grid for plotting
x, y, z = q.ParamGrid(rez=64)

# this can be used to view like this:
#from enthought.mayavi import mlab
#s = mlab.mesh(x, y, z, scalars=z**2)
#mlab.show()

x2 = x.ravel()
y2 = y.ravel()
z2 = z.ravel()
PlotPts = np.vstack((x2,y2,z2)).transpose()


j = maxlev - 0
TreeGrid.findStencilSets(j)
TreeGrid.buildListsFromStencilSets(j)
TreeGrid.findStencilsOnStencilSets(j, bdyfcn)
dx = TreeGrid.Levolve[j][0].dx
Lev = TreeGrid.Levolve[j]
Lex = TreeGrid.Lextend[j]
Grid = TreeGrid.Grids[j]
D = cpOps.buildDiffMatrix(Lev, Lex)
#D3xb,D3xf,D3yb,D3yf = g.buildDiffMatrixTempDxDyTest(j)
E = cpOps.buildExtensionMatrix(Lev, Lex)
Eplot = cpOps.buildEPlotMatrix(Grid, Lev, Lex, PlotPts, interp_degree=3)


# the implicit closest point method matrix.  Almost a product of D and
# E (see Macdonald and Ruuth 2009)
now = time()
M = cpOps.LinearDiagonalSplitting(D, E)
print "splitting took %.3g s" % (time()-now)


# keep a copy of the original extended precision matrix
#M_f96 = M.copy()
#M = M.astype(numpy.float64)


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
