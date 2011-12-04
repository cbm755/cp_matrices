"""
Closest point representation of a circular arc.

A circular arc in the x--y plane

Part of a 2D circle with center, radius as specified.  th1 and th2
specify the start and stop angles.

TODO: how to deal with pi for longdouble/float96?
"""
from Surface import ShapeWithBdy

from coordinate_transform import cart2pol, pol2cart
import numpy
from numpy import array as a
from numpy.linalg import norm
from numpy import pi

class CircularArc(ShapeWithBdy):
    def __init__(self, center=a([0.0, 0.0]), radius=1.0, angle1=pi/6, angle2=2*pi/3):
        self.center = center
        self.radius = radius
        self._dim = 2
        # TODO: this bounding box isn't great
        self._bb = [center - radius, center + radius]
        self._hasParam = True
        self.angle1 = angle1
        self.angle2 = angle2

    def closestPointToCartesian(self, xx):
        x,y = xx - self.center
        th,r = cart2pol(x, y)

        if ((th >= self.angle1) and (th <= self.angle2)):
            (cpx,cpy) = pol2cart(th, self.radius)
            bdy = 0
            dist = norm( a([x, y]) - a([cpx,cpy]), 2)
        else:
            # check the two end points
            (cpx,cpy) = pol2cart(self.angle1, self.radius)
            dist = norm( a([x, y]) - a([cpx,cpy]), 2)
            bdy = 1
            (cpx2,cpy2) = pol2cart(self.angle2, self.radius)
            dist2 = norm( a([x, y]) - a([cpx2,cpy2]), 2)
            if (dist2 < dist):
                cpx = cpx2
                cpy = cpy2
                dist = dist2
                bdy = 2

        cp = self.center + a([cpx,cpy])
        return cp, dist, bdy, {}


    cp = closestPointToCartesian

    def ParamGrid(self, rez=256):
        """
        Parameritized form (for plotting)
        """
        from numpy import linspace,pi,real,imag,exp

        th = linspace(self.angle1, self.angle2, num=rez, endpoint=True)
        circ = self.radius * numpy.exp(1j*th)
        X = real(circ) + self.center[0]
        Y = imag(circ) + self.center[1]
        #plot(real(circ), imag(circ), 'k-');
        #XtraPts = numpy.vstack((real(circ),imag(circ))).transpose()
        return X,Y
