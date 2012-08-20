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
        if self._dim == 2 or self._dim == 3:
            self._hasParam = True
        #self._dim = self._center.shape[0]
        # TODO: superclass knows about dimension?
        #super(Hemisphere, self).__init__()

    def closestPointToCartesianVec(self, x):
        x = np.atleast_2d(x)
        r = np.sqrt(np.sum((x - self._center)**2, axis=1))
        xm = x.copy()
        r_eq_0 = r == 0
        xm[r_eq_0, 0] = self._center[0] + 1
        r[r_eq_0] = 1
        dim = self._dim
        cpx = np.zeros(x.shape)
        bdy = np.zeros(x.shape[0])
        below_semicircle = xm[:, -1] < self._center[-1]
        if dim == 2:
            left = (xm[:, 0] <= self._center[0]) & below_semicircle
            cpx[left] = self._center - np.array([self._radius, 0.])
            bdy[left] = 1
            right = ~left & below_semicircle
            cpx[right] = self._center + np.array([self._radius, 0.])
            bdy[right] = 2
        elif dim == 3:
            # Here it doesn't matter if xm[below_semicircle, *] has
            # shape (0,) because adding a scalar "works" (ie, xx has
            # also shape (0,))
            xx = xm[below_semicircle, 0] - self._center[0]
            yy = xm[below_semicircle, 1] - self._center[1]
            th, _ = cart2pol(xx, yy)
            xx, yy = pol2cart(th, self._radius)
            cpx[below_semicircle] = (np.column_stack((xx,
                                                      yy ,
                                                      np.zeros(xx.shape))) +
                                     self._center)
            bdy[below_semicircle] = 1
        else:
            raise NotImplementedError(
                'Dim {0} and higher not implemented'.format(dim)
                )
        rest = ~below_semicircle
        
        # This is needed, if xm[rest] is empty (ie shape (0,)) trying
        # to add an array results in a ValueError: operands could not
        # be broadcast together
        if rest.any():
            c = (self._radius / r[rest])[:, np.newaxis]
            cpx[rest] = c * (xm[rest] - self._center) + self._center
            bdy[rest] = 0
        dist = np.sqrt(np.sum((x - cpx)**2, axis=1))
        return cpx, dist, bdy, {}
        
    def closestPointToCartesian(self, x):
        r = np.linalg.norm(x - self._center, 2)
        if r==0:
            tx = self._center[0] + 1.0
            r = 1.0
        else:
            tx = x[0]
        xm = x.copy()
        xm[0] = tx

        dim = self._dim

        if xm[-1] < self._center[-1]:
            # "below" semicircle,
            if dim == 2:
                if xm[0] <= self._center[0]:  # left
                    cpx = self._center - np.array([self._radius, 0.0])
                    bdy = 1
                else:  # right
                    cpx = self._center + np.array([self._radius, 0.0])
                    bdy = 2
            elif dim == 3:
                xx = xm[0] - self._center[0]
                yy = xm[1] - self._center[1]
                th, r = cart2pol(xx,yy)
                xx, yy = pol2cart(th, self._radius)
                cpx = self._center + np.array([xx,yy,0])
                bdy = 1
            else:
                # general case possible?
                raise NotImplementedError(
                    'Dim {0} and higher not implemented'.format(dim)
                    )
        else:
            # at or "above" semicircle: same as for whole circle
            c = self._radius / r
            cpx = c*(xm-self._center) + self._center
            bdy = 0

        dist = np.linalg.norm(cpx - x, 2)
        return cpx, dist, bdy, {}

    def closestPointToCartesianGrid(self, grid):
        """ Hack to accept grid arrays, since np.vectorize doesn't
        work for class methods (see
        mail.scipy.org./pipermail/numpy-discussion/2007-February/026060.html

        FIXME It looses the actual closestpoint, and bdy"""
        ret = np.empty(grid[0].shape).ravel()
        gr = grid.reshape((3, -1))
        for i in xrange(len(ret)):
            ret[i] = self.closestPointToCartesian(gr[:, i])[1]
        return None, ret.reshape(grid[0].shape), None, None
        
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
        return x + cen[0], y + cen[1], z + cen[2]


class Semicircle(Hemisphere):
    """
    Closest point represenation of a semicircle.
    """
    def __init__(self, center=np.array([0.0, 0.0]), radius=1.0):
        Hemisphere.__init__(self, center=center, radius=radius)

    def ParamGrid(self, rez=256):
        """
        Parameterized form (for plotting).
        """
        th = np.linspace(0, np.pi, rez)
        circ = self._radius * np.exp(1j*th)
        X = np.real(circ) + self._center[0]
        Y = np.imag(circ) + self._center[1]
        return X, Y
