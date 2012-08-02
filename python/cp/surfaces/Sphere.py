"""
Closest point function for a sphere.

Works in multidimensions.
"""
from Surface import Surface

from numpy import array as a
from scipy.linalg import norm
from scipy import sum as spsum
from scipy import sqrt,where,zeros


class Sphere(Surface):
    def __init__(self, center=a([0.0, 0.0, 0.0]), radius=1.0):
        self._center = center
        self._radius = radius
        self._dim = len(center)
        self._bb = [center - radius, center + radius]
        if ((self._dim == 2) or (self._dim == 3)):
            self._hasParam = True
        # TODO: superclass knows about dimension?
        #super(Hemisphere, self).__init__()

    def closestPointToCartesian(self, x):
        """
        This implementation works in general dimensions.  But not very
        "clean" and one has to be careful with the details (copy is
        necessary for the computation of distance for example)
        """
        SURF_CEN = self._center
        SURF_SR = self._radius
        
        if x.ndim == 1:
            cpx = x-SURF_CEN
            r = norm(cpx)
            if r == 0:
                cpx[0] += 1
                r = 1
            cpx = SURF_SR/r*cpx + SURF_CEN
            d = norm(cpx-x)
            return cpx,d,0,{}
            
        cpx = x-SURF_CEN
        r = sqrt(spsum(pow(cpx,2),axis=1))
        if (r.any()==0):
            ind = where(r == 0)
            cpx[ind,0] = 1
            r[ind] = 1


        r = SURF_SR / r
        for dim in range(x.shape[1]):
            cpx[:,dim] *= r
        cpx += SURF_CEN
        dist = sqrt(spsum(pow(cpx - x, 2),axis=1))
        return cpx, dist, zeros(x.ndim), {}

    cp = closestPointToCartesian

    def ParamGrid(self, rez=20):
        """ Return a mesh, for example for plotting with mlab
        """
        #import numpy as np
        from numpy import pi,sin,cos,outer,ones,size,linspace
        rad = self._radius
        cen = self._center
        # parametric variables
        #u=np.r_[0:2*pi:10j]
        #v=np.r_[0:0.5*pi:10j]
        th = linspace(0, 2*pi, num=2*rez, endpoint=True)
        phi = linspace(0, 1*pi, num=rez, endpoint=True)
        x = rad*outer(cos(th), sin(phi)) + cen[0]
        y = rad*outer(sin(th), sin(phi)) + cen[1]
        z = rad*outer(ones(size(th)),cos(phi)) + cen[2]
        # TODO: return a list?  In case we have multiple components
        return x,y,z
