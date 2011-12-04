# following the matlab approach of the cp_matrices code
#

#import numpy
import numpy as np
from numpy import array as a

# the reload here lets you use this in ipython as:
#   "run -i newapproach.py"
# even if changing the surfaces module.
import surfaces
reload(surfaces)

# there was already an old cpGrid object...
import cpGrid_new as cpGrid
reload(cpGrid)

import cp_ops_new as cp_ops
reload(cp_ops)

x1d,dx = np.linspace(-2,2,41,retstep=True)
y1d = x1d.copy()   # otherwise, its a pointer

xxg,yyg = np.meshgrid(x1d, y1d)

# vectors instead of 2D arrays
xx = xxg.flatten()
yy = yyg.flatten()
# TODO: make this a n x 2 matrix...
xy = np.hstack( (xx.reshape(xx.shape[0],1), yy.reshape(yy.shape[0],1)) )


#from closestpoint import coordinate_transform as tf
q = surfaces.Circle()

# wrap the object in CPBar, for accurately imposing boundary
# conditions on shapes with boundaries (like a hemisphere)
#q = closestpoint.CPBar(parent=q)

cpfun = q.closestPointVectorized

#th, r = tf.cart2pol(xx, yy)
#cpx, cpy = tf.pol2cart(th, 1.1)
#cpx, cpy = tf.pol2cart(th, self._radius)

#A = q.cp(xx,yy)
cpx,cpy,dist,bdy,xtra = cpfun(xx,yy)


# Banding: do calculation in a narrow band around the circle
dim = 2;  # dimension
p = 3;    # interpolation degree
# "band" is a vector of the indices of the points in the computation
# band.  The formula for bw is found in [Ruuth & Merriman 2008] and
# the 1.0001 is a safety factor.
bw = 1.0001*np.sqrt((dim-1)*((p+1)/2)**2 + ((1+(p+1)/2)**2))

# this comma is signficant
band, = np.nonzero(np.abs(dist) <= bw*dx)

cpx = cpx[band]
cpy = cpy[band]
dist = dist[band]
#bdy = bdy[band]
x = xx[band]
y = yy[band]
xy = xy[band,:]

g1 = cpGrid.cpGrid(x1d, y1d, dx)
g1.cpx = cpx
g1.cpy = cpy
g1.band = band
g1.dist = dist
g1.x = x
g1.y = y
g1.xy = xy

#IJ = cpGrid.i2s(xxg.shape,band)

# a grid nows how to convert subscripts to linear indices and vice-versa
g1.ij = g1.ind2sub(band)

# double check
#band_check = g1.sub2ind(g1.ij)
#print g1.band - band_check

#q.viztest()

# two codes for the cartesian laplacian
L  = cp_ops.buildDiffMatrix(g1)
L2  = cp_ops.buildDiffMatrixFast(g1)
L - L2

# make it square
LL = L[:,band]

cpxy = np.hstack( (cpx.reshape(cpx.shape[0],1), cpy.reshape(cpy.shape[0],1)) )


# also returns an inner band for a dual-banded code (not used here)
(E,band2) = cp_ops.buildExtensionMatrix(g1, cpxy, degreep=3)
