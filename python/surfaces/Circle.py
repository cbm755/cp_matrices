"""
Closest point representation of a circle.

A special case of Sphere but implemented differently because this
version is clearer to understand.
"""
from Surface import Surface

from coordinate_transform import cart2pol, pol2cart
import numpy
from numpy import array as a
from numpy.linalg import norm

class Circle(Surface):
    def __init__(self, center=a([0.0, 0.0]), radius=1.0):
        # TODO: could make it subclass of sphere and just call:
        # self.super(self, center=center, radius=radius)
        self._center = center
        self._radius = radius
        self._dim = 2
        self._bb = [center - radius, center + radius]
        self._hasParam = True

    def closestPointVectorized(self, x, y):
        th, r = cart2pol(x - self._center[0], y - self._center[1])
        cpx, cpy = pol2cart(th, self._radius)

        #dist = norm(xx - cp, 2)
        #dist = sqrt( (x-cpx)**2 + (y-cpy)**2 )
        sdist = r - self._radius
	cpx = cpx + self._center[0]
	cpy = cpy + self._center[1]
        return cpx, cpy, sdist, 0, {}

    def closestPointToCartesian(self, xx):
        # TODO: could probably be vectorized
        x,y = xx - self._center

        th, r = cart2pol(x, y)
        x, y = pol2cart(th, self._radius)
        cp = self._center + a([x,y])

        dist = norm(xx - cp, 2)
        return cp, dist, 0, {}

    cp = closestPointToCartesian

    def ParamGrid(self, rez=256):
        """
        Parameritized form (for plotting)
        """
        from numpy import linspace,pi,real,imag,exp

        th = linspace(0, 2*pi, num=rez, endpoint=True)
        circ = self._radius * numpy.exp(1j*th)
        X = real(circ) + self._center[0]
        Y = imag(circ) + self._center[1]
        #plot(real(circ), imag(circ), 'k-');
        #XtraPts = numpy.vstack((real(circ),imag(circ))).transpose()
        return X,Y
