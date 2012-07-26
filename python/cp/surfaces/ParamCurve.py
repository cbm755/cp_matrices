"""
Closest Point Representation of a parameterized curve, found by
optimization.

TODO: currently hardcoded to (x,y) = (t,cos(t)) in 2D, easy to
generalize this to 3D, (x,y,z) = (f(t),g(t),h(t))

Newton's method on derivative of distance squared works better than
straight optimization on distance squared.
"""
import numpy as np
from scipy.optimize import fminbound

from Surface import ShapeWithBdy


class ParamCurve(ShapeWithBdy):
    def __init__(self, end1=1.0/4, end2=4, f=None):
        l = end1
        r = end2
        self._l = end1
        self._r = end2
        self._dim = 2
        if f is None:
            def f(t):
                return np.cos(t)
        def mf(t):
            return -f(t)
        self._f = f

        self._tol = 1e-14
        tol = self._tol

        # Use optimization to find the bounding box
        (topt, fmin, ierr, numfunc) = fminbound(f,  l, r, xtol=tol,full_output=True,disp=1)
        (topt, fval, ierr, numfunc) = fminbound(mf,l , r, xtol=tol,full_output=True,disp=1)
        fmax = -fval
        self._bb = [np.array([l,fmin]), np.array([r,fmax]) ]

        self._hasParam = True


    def closestPointToCartesianOld(self, xx):
        """ this version doesn't detect endpts """
        x, y = xx
        f = self._f

        def d2(s,x,y):
            return (s-x)**2 + (f(s)-y)**2

        endpt1 = self._l
        endpt2 = self._r
        (sopt, fval, ierr, numfunc) = fminbound(d2, endpt1, endpt2, args=(x,y), \
                                                    xtol=self._tol,full_output=True,disp=3)
        cp = np.array([sopt,cos(sopt)])

        dist = np.sqrt(fval)
        dist2 = np.linalg.norm(xx-cp, 2)
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
            raise ValueError('Max iter exceeded')
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

        cp = np.array([tmin, self._f(tmin)])
        #dist2 = sqrt( (xx[0]-cp[0])**2 + (xx[1]-cp[1])**2 )
        dist = np.linalg.norm(xx - cp, 2)

        others = dict(param=tmin)
        return cp, dist, bdy, others


    def closestPointToCartesianColin(self, xx, verbose=0):
        x,y = xx
        f = self._f
        fp = lambda s: -np.sin(s)
        fpp = lambda s: -np.cos(s)
        endpt1 = self._l
        endpt2 = self._r

        def d2(s, x, y):
            return (s - x)**2 + (f(s) - y)**2
        #def d2(s, yy):
        #    return (s - yy[0])**2 + (f(s) - yy[1])**2
        # dist squared deriv
        def g(s, x, y):
            return 2*(s - x) + 2*(f(s) - y)*fp(s)
            #return 2*s - 2*x - 2*(cos(s) - y)*sin(s)
        def gp(s, x, y):
            return 2 + 2*(f(s) - y)*fpp(s) + 2*fp(s)*fp(s)
            #return 2 + 2*(sin(s))**2 - 2*(cos(s) - y)*cos(s)

        tmin = None
        ddmin = 1e42

        disp = 1
        if  (verbose > 12): disp=3
        t,dd,ierr,numfunc = fminbound(d2, endpt1, endpt2, args=(xx[0],xx[1]), \
                                          full_output=True, \
                                          xtol=self._tol, maxfun=5000, disp=disp)
        if ierr == 1:
            raise ValueError('max iter exceeded')
        tmin = t[0]  # for same reason t is a length 1 array
        ddmin = dd

        bdy = 0

        if ((abs(tmin - endpt1) < 1e-4) or (abs(tmin - endpt2) < 1e-4)):
            #print 'close to endpt, careful', x, y
            CloseToEndpt = True
        else:
            CloseToEndpt = False

        # tighten with newton's method:
        #print "tighten with newton's method", x, y

        s = tmin.astype('float96')
        #s = tmin
        n = 1
        while True:
            #print type(g(s,x,y))
            n = n + 1
            snew = s - g(s,x,y) / gp(s,x,y)
            if (verbose > 10):
                print s-snew
            if (abs(s - snew) < 5e-19):
                s = snew
                fail = False
                break
            if (n > 1000):
                fail = True
                break
            s = snew
        if fail:
            if not CloseToEndpt:
                print 'too many iterations, not close to endpoint'
                print x,y
                print tmin, s, tmin-s



        if CloseToEndpt:
            dd2 = d2(s,x,y)
            if ((s > endpt1) and (s < endpt2) and (dd2 <= dd)):
                print "newton says not endpoint"
                print "(x,y)=", x, y, "(tmin,s,diff)=", tmin, s, tmin - s
                print "dd=", dd, "dd2=", d2(s,x,y)
        else:
            if (abs(tmin - s) > 1e-7):
                print "warning"
                print "CloseToEndpt=", CloseToEndpt
                print "(x,y)=", x, y, "(tmin,s,diff)=", tmin, s, tmin - s
            if ((s < endpt1) or (s > endpt2)):
                print "** warning endpt **"
                print x,y
                print tmin, s, tmin-s
            
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

        cp = np.array([tmin, self._f(tmin)])
        #dist2 = sqrt( (xx[0]-cp[0])**2 + (xx[1]-cp[1])**2 )
        dist = np.linalg.norm(xx - cp, 2)

        others = dict(param=tmin)
        return cp, dist, bdy, others



    def closestPointToCartesian_Newton(self, xx, verbose=0):
        """ Newton's method """
        x, y = xx
        f = self._f
        fp = lambda s: -np.sin(s)
        fpp = lambda s: -np.cos(s)
        endpt1 = self._l
        endpt2 = self._r

        def d2(s, x, y):
            return (s - x)**2 + (f(s) - y)**2
        #def d2(s, yy):
        #    return (s - yy[0])**2 + (f(s) - yy[1])**2
        # dist squared deriv
        def g(s, x, y):
            return 2*(s - x) + 2*(f(s) - y)*fp(s)
            #return 2*s - 2*x - 2*(cos(s) - y)*sin(s)
        def gp(s, x, y):
            return 2 + 2*(f(s) - y)*fpp(s) + 2*fp(s)*fp(s)
            #return 2 + 2*(sin(s))**2 - 2*(cos(s) - y)*cos(s)

        # Todo: time it, balance this with the newton solve
        ss = np.linspace(endpt1, endpt2, 1000)
        dd = d2(ss, x, y)
        t_guess = ss[dd.argmin()]
        dd_guess = dd[dd.argmin()]
        t_guess= t_guess.astype('float96')


        # tighten with newton's method:
        #print "tighten with newton's method", x, y


        s = t_guess
        n = 1
        outsideCounter = 0
        while True:
            #print type(g(s,x,y))
            n = n + 1
            snew = s - g(s,x,y) / gp(s,x,y)
            # this escape tolerance should be set somewhere
            if (snew - endpt1 < -10):
                outsideCounter += 1
            if (snew - endpt2 > 10):
                outsideCounter += 1
            if (verbose > 10):
                print s-snew
            if (abs(s - snew) < 4e-19):  # 4*macheps, here using extended precision, TODO
                fail = False
                #print "converged: (x,y,s,n)=",x,y,snew,n
                #print n
                #print s, snew
                break
            if (outsideCounter >= 10):
                fail = True
                print "**** probably diverging: (t_guess,n,s,snew)=", t_guess,n,s,snew, "endpts=", endpt1, endpt2
                break
            if (n > 255):
                #fail = True
                #break
                #raise NameError('too many iterations')
                #print 'too many iterations', n, snew, x, y
                fail = True
                print "XXXX Too many iterations: (t_guess,n,s,snew)=", t_guess,n,s,snew, "endpts=", endpt1, endpt2
                break
            s = snew

        #if (fail):
        #    if not CloseToEndpt:
        #        print 'too many iterations, not close to endpoint'
        #        print x,y
        #        print tmin, s, tmin-s

        #if fail:
        #    ddmin = 10000
        #    tmin = -1000
        #else:
        #    ddmin = d2(s, x, y)
        #    tmin = snew

        if fail:
            #dd_guess = dd[dd.argmin()]
            dd1 = d2(endpt1, x, y)
            dd2 = d2(endpt2, x, y)
            if ((dd_guess < dd1 - 1e-14) and (dd_guess < dd2 - 1e-14)):
                print dd_guess, dd1, dd2
                print dd_guess - dd1, dd_guess - dd2
                raise NameError('really?')
            if (dd1 <= dd2):
                bdy = 1
                ddmin = dd1
                tmin = endpt1
            else:
                bdy = 2
                ddmin = dd2
                tmin = endpt2
        else:
            # todo: one failure case leads to better end points detection
            if snew <= endpt1:
                bdy = 1
                tmin = endpt1
                ddmin = d2(tmin, x, y)
            elif snew >= endpt2:
                bdy = 2
                tmin = endpt2
                ddmin = d2(tmin, x, y)
            else:
                bdy = 0
                tmin = snew
                ddmin = d2(tmin, x, y)

        # TODO: compute the old way and check
        if (1==0):
            tmin = None
            ddmin = 1e42

            disp = 1
            if  (verbose > 12): disp=3
            t,dd,ierr,numfunc = fminbound(d2, endpt1, endpt2, args=(xx[0],xx[1]), \
                                          full_output=True, \
                                          xtol=self._tol, maxfun=5000, disp=disp)
            if ierr == 1:
                raise NameError('max iter exceeded')
            tmin = t[0]  # for same reason t is a length 1 array
            ddmin = dd

            bdy = 0

            if ((abs(tmin - endpt1) < 1e-4) or (abs(tmin - endpt2) < 1e-4)):
                #print 'close to endpt, careful', x, y
                CloseToEndpt = True
            else:
                CloseToEndpt = False


        cp = np.array([tmin, self._f(tmin)])
        #dist2 = sqrt( (xx[0]-cp[0])**2 + (xx[1]-cp[1])**2 )
        dist = np.linalg.norm(xx - cp, 2)

        others = dict(param=tmin)
        return cp, dist, bdy, others

    cp = closestPointToCartesian_Newton



    def ParamGrid(self, rez=256):
        """
        Parameritized form (for plotting)
        """
        th = np.linspace(self._l, self._r, rez)
        X = th
        Y = np.cos(th)
        return X,Y
