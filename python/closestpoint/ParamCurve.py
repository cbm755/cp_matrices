"""
Closest Point Representation of a parameterized curve, found by
optimization.

TODO: currently hardcoded to (x,y) = (t,cos(t)) in 2D, easy to
generalize this to 3D, (x,y,z) = (f(t),g(t),h(t))

TODO: tolerances?  My tests suggest putting them to 1e-15 doesn't
improve the result.
"""
from ClosestPoint import ShapeWithBdy

from numpy import array as a
from numpy import sqrt, linspace, cos, pi
from numpy.linalg import norm
from scipy.optimize import fminbound

class ParamCurve(ShapeWithBdy):
    def __init__(self, end1=0, end2=3*pi/2, f=None):
        l = end1
        r = end2
        self._l = end1
        self._r = end2
        self._dim = 2
        if f==None:
            def f(t):
                return cos(t)
        def mf(t):
            return -f(t)
        self._f = f

        self._tol = 1e-12
        tol = self._tol

        # Use optimization to find the bounding box
        (topt, fmin, ierr, numfunc) = fminbound(f,  l, r, xtol=tol,full_output=True,disp=1)
        (topt, fval, ierr, numfunc) = fminbound(mf,l , r, xtol=tol,full_output=True,disp=1)
        fmax = -fval
        self._bb = [ a([l,fmin]), a([r,fmax]) ]

        self._hasParam = True


    def closestPointToCartesianOld(self, xx):
        x,y = xx
        f = self._f

        def d22(s,x,y):
            return (s-x)**2 + (f(s)-y)**2

        endpt1 = self._l
        endpt2 = self._r
        (sopt, fval, ierr, numfunc) = fminbound(d22, endpt1, endpt2, args=(x,y), \
                                                    xtol=1e-12,full_output=True,disp=3)
        cp = a([sopt,cos(sopt)])

        dist = sqrt(fval)
        dist2 = norm(xx-cp, 2)
        print dist-dist2

        return cp, dist2, 0


    def closestPointToCartesian(self, xx):
        #x,y = xx
        f = self._f
        endpt1 = self._l
        endpt2 = self._r

        print 'hello'

        def d2(s, x, y):
            return (s - x)**2 + (f(s) - y)**2

        tmin = None
        ddmin = 1e42

        t,dd,ierr,numfunc = fminbound(d2, endpt1, endpt2, args=(xx), \
                                          full_output=True, \
                                          xtol=1e-12, maxfun=5000, disp=1)
        # TODO: some error checking?
        tmin = t[0]  # for same reason t is a length 1 array
        ddmin = dd

        bdy = 0

        dd = d2(endpt1, xx[0], xx[1])
        # TODO: add some small multiple of macheps?
        if dd <= ddmin:
            bdy = 1
            ddmin = dd
            tmin = endpt1

        dd = d2(endpt2, xx[0], xx[1])
        if dd <= ddmin:
            bdy = 2
            ddmin = dd
            tmin = endpt2

        #cp = a([0.0, 0.0])
        #(x,y) = scipy.interpolate.splev(tmin, sp.tck)
        # TODO: hardcoded for 2D

        cp = a([tmin, self._f(tmin)])
        dist = sqrt( (xx[0]-cp[0])**2 + (xx[1]-cp[1])**2 )
        dist2 = norm(xx - cp, 2)
        print dist-dist2, bdy, tmin

        # TODO: change all the codes to return an optional list
        #return cp, dist, bdy, (tmin)
        return cp, dist, bdy

    cp = closestPointToCartesian



    def ParamGrid(self, rez=256):
        """
        Parameritized form (for plotting)
        """
        th = linspace(self._l, self._r, num=rez, endpoint=True)
        X = th
        Y = cos(th)
        return X,Y

