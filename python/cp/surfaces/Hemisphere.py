"""
Closest point representation for hemisphere and semicircle.

Should work in multidimensions: the hemi part happens in the last dimension.
Not necessarily all implemented yet

TODO: WIP
"""
import numpy as np


from Surface import ShapeWithBdy
from coordinate_transform import cart2pol, pol2cart, sph2cart

class Hemisphere(ShapeWithBdy):
    def __init__(self, center=np.array([0.0, 0.0, 0.0]), radius=1.0):
        self._center = center
        self._radius = radius
        self._dim = len(center)
        ll = center - radius
        rr = center + radius
        ll[-1] = center[-1]
        self._bb = [ll, rr]
        if (self._dim == 2) or (self._dim == 3):
            self._hasParam = True
        #self._dim = self._center.shape[0]
        # TODO: superclass knows about dimension?
        #super(Hemisphere, self).__init__()

    def closestPointToCartesian(self, x):
        SURF_CEN = self._center
        SURF_SR = self._radius

        r = np.linalg.norm(x - SURF_CEN, 2)
        if r==0:
            tx = SURF_CEN[0] + 1.0
            r = 1.0
        else:
            tx = x[0]
        xm = x.copy()
        xm[0] = tx

        #dim = len(x)
        dim = self._dim

        if xm[-1] < SURF_CEN[-1]:
            # "below" semicircle,
            if dim == 2:
                if xm[0] <= SURF_CEN[0]:  # left
                    cpx = SURF_CEN - np.array([SURF_SR, 0.0])
                    bdy = 1
                else:  # right
                    cpx = SURF_CEN + np.array([SURF_SR, 0.0])
                    bdy = 2
            if dim == 3:
                xx = xm[0] - SURF_CEN[0]
                yy = xm[1] - SURF_CEN[1]
                th, r = cart2pol(xx,yy)
                xx, yy = pol2cart(th,SURF_SR)
                cpx = SURF_CEN + np.array([xx,yy,0])
                bdy = 1
            if dim >= 4:
                # general case possible?
                raise NotImplementedError('4D and higher not implemented')
        else:
            # at or "above" semicircle: same as for whole circle
            c = SURF_SR / r
            cpx = c*(xm-SURF_CEN) + SURF_CEN
            bdy = 0

        dist = np.linalg.norm(cpx - x, 2)
        return cpx, dist, bdy, {}

    cp = closestPointToCartesian


    def ParamGrid(self, rez=10):
        """
        Return a mesh, for example for plotting with mlab.
        """
        rad = self._radius
        cen = self._center
        # parametric variables
        #u=np.r_[0:2*pi:10j]
        #v=np.r_[0:0.5*pi:10j]
        th = np.linspace(0, 2*np.pi, 4*rez)
        phi = np.linspace(0, 0.5*np.pi, rez)
        TH, PHI = np.meshgrid(th, phi)
        x, y, z = sph2cart(TH, PHI, rad)
        # TODO: return a list?  In case we have multiple components
        return x.ravel() + cen[0], y.ravel() + cen[1], z.ravel() + cen[2]


class Semicircle(Hemisphere):
    """
    Closest point represenation of a semicircle.
    """
    def __init__(self, center=np.array([0.0, 0.0]), radius=1.0):
        #print super(self)
        #self.super(self, center=center, radius=radius)
        # call the superclass constructor
        Hemisphere.__init__(self, center=center, radius=radius)
        #self._hasParam = True

    def ParamGrid(self, rez=256):
        """
        Parameterized form (for plotting).
        """
        th = np.linspace(0, 1*np.pi, num=rez, endpoint=True)
        circ = self._radius * np.exp(1j*th)
        X = np.real(circ) + self._center[0]
        Y = np.imag(circ) + self._center[1]
        #plot(real(circ), imag(circ), 'k-');
        #XtraPts = numpy.vstack((real(circ),imag(circ))).transpose()
        return X, Y



#if __name__ == "__main__":
#    print "running as a script, do some tests"
#   test()
