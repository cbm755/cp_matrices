# following the matlab approach of the cp_matrices code
#

import numpy as np
import numpy

x1d,dx = np.linspace(-2,2,11,retstep=True)
y1d = x1d.copy()   # otherwise, its a pointer

xxg,yyg = np.meshgrid(x1d, y1d)

# vectors instead of 2D arrays
xx = xxg.flatten()
yy = yyg.flatten()
# or
xy = np.hstack((xx.T, yy.T))

import surfaces
reload(surfaces)
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

#q.viztest()

IJ = i2s(xxg.shape,band)



def i2s(sz,ind):
    """works like ind2sub but with v = s2i(sz,ind) instead of
    [v1,v2,v3,...] = ind2sub(sz,ind)
    works along dimension 2 (i.e. rows specify the index)
    CAUTION no error checking (if index out of range)"""

    ind = numpy.array(ind) + 1
    sz = numpy.array(sz)
    sub=numpy.zeros([ind.size,sz.size]).astype(numpy.int);

    for i in range(len(sz)-1):

        sub[:,i] = numpy.mod(ind,sz[i]);
        sub[sub[:,i]==0,i] = sz[i];
        ind = (ind - sub[:,i])/sz[i] + 1;

    sub[:,i+1] = ind;
    sub = sub - 1

    return sub



def s2i(sz,v):
    """ works like sub2ind but with s2i(sz,[v1,v2,v3,...]) instead of
    sub2ind(sz,v1,v2,v3,...)
    works along dimension 2 (i.e. rows specify the indices)
    CAUTION no error checking (if index out of range)"""

    index = numpy.array(v[:,0]).astype(numpy.int)
    sc=1
    for i in range(1,len(sz)):
        sc = sc*(sz[i-1]);
        index = index + sc*(v[:,i]);

    return index
