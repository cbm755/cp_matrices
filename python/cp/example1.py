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

import surfaces.coordinate_transform
from surfaces import coordinate_transform as ctf

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
# make this a n x 2 matrix...
xy = np.hstack( (xx.reshape(xx.shape[0],1), yy.reshape(yy.shape[0],1)) )


# a surface object, these know how to compute their own cp's
q = surfaces.Circle()

# wrap the object in CPBar, for accurately imposing boundary
# conditions on shapes with boundaries (like a hemisphere)
#q = closestpoint.CPBar(parent=q)

# a pointer-to-function
cpfun = q.closestPointVectorized

cpx,cpy,dist,bdy,xtra = cpfun(xx,yy)

# or can compute the cp data manually:
#th, r = ctf.cart2pol(xx, yy)
#cpx, cpy = ctf.pol2cart(th, 1.1)



# no banding in this code, so just set band to everything
band = range(0, len(x1d)*len(y1d))

#x = xx[band]
#y = yy[band]
#xy = xy[band,:]
x = xx
y = yy

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
LL  = cp_ops.buildDiffMatrix(g1)
LL2  = cp_ops.buildDiffMatrixFast(g1)
if (LL - LL2).nnz > 0:
    raise NameError('two laplacian codes produced different results')

# LL operators on the meshgrid, reduce to the band
L = LL[:,band]

cpxy = np.hstack( (cpx.reshape(cpx.shape[0],1), cpy.reshape(cpy.shape[0],1)) )

# also returns an inner band for a dual-banded code (not used here)
(EE,band2) = cp_ops.buildExtensionMatrix(g1, cpxy, degreep=3)

E = EE[:,band]

(th,r) = ctf.cart2pol(cpx,cpy)

u0 = np.cos(th)

Tf = 0.25
dt = 0.2*dx**2

nt = np.ceil(Tf/dt).astype(int)
dt = Tf / nt

u = u0
for k in range(0, nt):
    t = (k+1)*dt
    unew = u + dt*(L*u)
    u = E*unew;

uex = np.exp(-t) * np.cos(th)

err = (abs(u-uex)).max()
print err

import pylab

plot = pylab.plot


plot(th, u0, 'bx', label='IC')

plot(th, E*u, 'ro', label='u num')

plot(th, uex, 'gx', label='u exact')

pylab.xlabel(r"$\theta$")
pylab.ylabel(r"$u$")
pylab.legend()
#legend( ('label1', 'label2', 'label3'), loc='upper left')
# ipython --pylab
# then this is non-blocking
pylab.show()

