"""
Closest point function for a sphere.

Works in multidimensions.
"""
from ClosestPoint import ClosestPoint

from numpy import array as a
from scipy.linalg import norm

class Sphere(ClosestPoint):
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

        r = norm(x - SURF_CEN, 2)
        if (r==0):
            tx = SURF_CEN[0] + 1.0
            r = 1.0
        else:
            tx = x[0]
        xm = x.copy()
        xm[0] = tx

        c = SURF_SR / r
        cpx = c*(xm-SURF_CEN) + SURF_CEN
        dist = norm(cpx - x, 2)
        return cpx, dist, 0, {}

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
