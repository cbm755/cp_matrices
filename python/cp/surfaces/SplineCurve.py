"""
Closest Point Representation of a spline curve, found by
optimization using Newton's method.

TODO: this is periodic for now...

TODO: code needs some cleanup, there is still stuff left from a
optimization code with endpoints.
"""
from Surface import Surface

import numpy as np

from scipy.interpolate import splprep, splev
from scipy.optimize import fminbound
import scipy.interpolate as si
import scipy.integrate

class SplineCurve(Surface):
    def __init__(self, sc=1.68882661023481):  #, preset='egg'):
        #scale_matlab = 1.6917222947492

        # this value is determined experimentally to make it 2*pi
        # arclength (seems to differ from matlab (different
        # smoothing?)
        scale=1.6888266102348

        pts = scale*np.array([
                             [0, 0.7],
                             [0.4, 0.1],
                             [0.2, -0.7],
                             [-0.4, -0.5],
                             [-0.4, 0.2],
                             [0, 0.7]])
        x = pts[:,0]
        y = pts[:,1]
        pts2 = [x, y]

        (tck, u) = splprep(pts2, s=0, per=True)

        # find the arclength
        def integrand(s):
            return np.sqrt(np.sum(np.array(si.splev(s, tck, der=1))**2  ))

        arclen,errest = scipy.integrate.quad(integrand,0,1,epsabs=1e-15, epsrel=5e-14, limit=500)

        print "arclen,errst,err", arclen, errest, 2*np.pi-arclen

        self.tck = tck
        self.u = u
        self.paramBounds = [u[0], u[-1]]
        self._dim = 2

        # find an approximate bounding box
        s = np.linspace(u[0], u[-1], 1000)
        x,y = splev(s, self.tck)
        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()
        width = xmax-xmin
        height = ymax-ymin
        # add/subtract 1% of the width for a little padding
        self._bb = [np.array([xmin-width/100,ymin-height/100]),
                    np.array([xmax+width/100,ymax+height/100])]

        self._hasParam = True


    def closestPointToCartesian(self, xx, verbose=0):
        """Using Newton's method """
        x, y = xx
        endpt1 = self.paramBounds[0]
        endpt2 = self.paramBounds[1]

        def d2(s, X, Y):
            x, y = splev(s, self.tck)
            return (x - X)**2 + (y - Y)**2
        # dist squared deriv
        def g(s, X, Y):
            x, y = splev(s, self.tck)
            xp, yp = splev(s, self.tck, der=1)
            return 2*(x - X)*xp + 2*(y - Y)*yp
        def gp(s, X, Y):
            x, y = splev(s, self.tck)
            xp, yp = splev(s, self.tck, der=1)
            xpp, ypp = splev(s, self.tck, der=2)
            return 2*(xp**2) + 2*(x-X)*(xpp)  +  2*(yp)**2 + 2*(y-Y)*(ypp)

        # Todo: time it, balance this with the newton solve
        ss = np.linspace(endpt1, endpt2, 1000)
        dd = d2(ss, x, y)
        t_guess = ss[dd.argmin()]
        dd_guess = dd[dd.argmin()]
        # TODO: splev won't work with float96?
        #t_guess= t_guess.astype('float96')


        def argPeriodic(v, P=2*np.pi):
            """ Return the principal value the argument. """
            A = v/P
            B = A - np.floor(A)
            return B*P

        # tighten with newton's method:
        #print "tighten with newton's method", x, y
        s = t_guess
        n = 1
        outsideCounter = 0
        while True:
            #print type(g(s,x,y))
            n = n + 1
            snew = s - g(s,x,y) / gp(s,x,y)
            # careful of going outside the periodic boundaries
            # TODO: bug here: assumes paramBounds[0] == 0, true for egg
            snew = argPeriodic(snew, P=self.paramBounds[1])

            # this escape tolerance should be set somewhere
            if (snew - endpt1 < -10):
                outsideCounter += 1
            if (snew - endpt2 > 10):
                outsideCounter += 1
            if (verbose > 10):
                print s-snew
            if (abs(s - snew) < 5e-16):  # TODO: hardcoded tolerance here
                fail = False
                #print "converged: (x,y,s,n)=",x,y,snew,n
                #print n
                #print s, snew
                break
            if (outsideCounter >= 10):
                fail = True
                print "**** probably diverging: (t_guess,n,s,snew)=", t_guess,n,s,snew, "endpts=", endpt1, endpt2
                break
            if (n > 1000):
                #fail = True
                #break
                #raise NameError('too many iterations')
                #print 'too many iterations', n, snew, x, y
                fail = True
                print "XXXX Too many iterations: (t_guess,n,s,snew,diff)=", t_guess,n,s,snew,s-snew, "endpts=", endpt1, endpt2
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
            t,dd,ierr,numfunc = fminbound(d2, endpt1, endpt2, args=(xx[0],xx[1]),
                                          full_output=True,
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

        cp = np.array(splev(tmin, self.tck))
        dist = np.linalg.norm(xx - cp, 2)

        others = dict(param=tmin)
        return cp, dist, bdy, others

    cp = closestPointToCartesian



    def ParamGrid(self, rez=256):
        """
        Parameritized form (for plotting)
        """
        s = np.linspace(self.paramBounds[0], self.paramBounds[1], rez)
        x,y = splev(s, self.tck)
        return x, y
