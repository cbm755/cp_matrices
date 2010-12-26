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
from numpy import sqrt, linspace, pi  #, cos
from numpy import cos as npcos
from math import cos

from numpy.linalg import norm
from scipy.optimize import fminbound

class ParamCurve(ShapeWithBdy):
    def __init__(self, end1=1.0/4, end2=4, f=None):
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

        self._tol = 1e-14
        tol = self._tol

        # Use optimization to find the bounding box
        (topt, fmin, ierr, numfunc) = fminbound(f,  l, r, xtol=tol,full_output=True,disp=1)
        (topt, fval, ierr, numfunc) = fminbound(mf,l , r, xtol=tol,full_output=True,disp=1)
        fmax = -fval
        self._bb = [ a([l,fmin]), a([r,fmax]) ]

        self._hasParam = True


    def closestPointToCartesianOld(self, xx):
        """ this version doesn't detect endpts """
        x,y = xx
        f = self._f

        def d22(s,x,y):
            return (s-x)**2 + (f(s)-y)**2

        endpt1 = self._l
        endpt2 = self._r
        (sopt, fval, ierr, numfunc) = fminbound(d22, endpt1, endpt2, args=(x,y), \
                                                    xtol=self._tol,full_output=True,disp=3)
        cp = a([sopt,cos(sopt)])

        dist = sqrt(fval)
        dist2 = norm(xx-cp, 2)
        print dist-dist2

        others = dict(param=sopt)
        return cp, dist2, 0, others


    def closestPointToCartesian(self, xx):
        #x,y = xx
        f = self._f
        endpt1 = self._l
        endpt2 = self._r

        def d2(s, x, y):
            return (s - x)**2 + (f(s) - y)**2
        #def d2(s, yy):
        #    return (s - yy[0])**2 + (f(s) - yy[1])**2

        tmin = None
        ddmin = 1e42

        t,dd,ierr,numfunc = fminbound(d2, endpt1, endpt2, args=(xx[0],xx[1]), \
                                          full_output=True, \
                                          xtol=self._tol, maxfun=5000, disp=1)
        if ierr == 1:
            raise nameError('max iter exceeded')
        tmin = t[0]  # for same reason t is a length 1 array
        ddmin = dd

        bdy = 0

        dd = d2(endpt1, xx[0], xx[1])
        # TODO: add some small multiple of macheps?
        if dd <= (ddmin + 1e-10):
            bdy = 1
            ddmin = dd
            tmin = endpt1

        dd = d2(endpt2, xx[0], xx[1])
        if dd <= (ddmin + 1e-10):
            bdy = 2
            ddmin = dd
            tmin = endpt2

        #cp = a([0.0, 0.0])
        #(x,y) = scipy.interpolate.splev(tmin, sp.tck)
        # TODO: hardcoded for 2D

        cp = a([tmin, self._f(tmin)])
        #dist2 = sqrt( (xx[0]-cp[0])**2 + (xx[1]-cp[1])**2 )
        dist = norm(xx - cp, 2)

        others = dict(param=tmin)
        return cp, dist, bdy, others

    cp = closestPointToCartesian



    def ParamGrid(self, rez=256):
        """
        Parameritized form (for plotting)
        """
        th = linspace(self._l, self._r, num=rez, endpoint=True)
        X = th
        Y = npcos(th)
        return X,Y

