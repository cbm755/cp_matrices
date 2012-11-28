"""
Closest Point Representation of a parameterized closed curve in 2D, found by
Newton's method.


xs, ys : parameterized curve
xp, yp : first derivative of curve function
xpp, ypp : second derivative of curve function
endpt1, endpt2 : endpts for curve with boundaries / endvalues of parameter

"""

import numpy as np
import warnings
from scipy.optimize import fminbound
from Surface import Surface


class ParamCurveClosed(Surface):
    def __init__(self, xs, ys, xp, yp, xpp, ypp, endpt1=0.0, endpt2=1.0):
        
        self.dim = 2
        self._xs = xs
        self._ys = ys
        self._xp = xp
        self._yp = yp
        self._xpp = xpp
        self._ypp = ypp
        self._l = endpt1
        self._r = endpt2
 
        self.bounding_box = self.find_bounding_box()
        self._hasParam = True


    def find_bounding_box(self):
        """
        Use optimization to find the bounding box
        """
        l = self._l
        r = self._r
        xs = self._xs
        ys = self._ys
        tol = 1e-14
        def mxs(t):
            return -xs(t)
        def mys(t):
            return -ys(t)
        (sopt, xmin, ierr, numfunc) = fminbound(xs, l, r, xtol=tol,full_output=True,disp=1)
        (sopt, xval, ierr, numfunc) = fminbound(mxs, l , r, xtol=tol,full_output=True,disp=1)
        (sopt, ymin, ierr, numfunc) = fminbound(ys, l, r, xtol=tol,full_output=True,disp=1)
        (sopt, yval, ierr, numfunc) = fminbound(mys, l , r, xtol=tol,full_output=True,disp=1)
        xmax = -xval
        ymax = -yval
        return [np.array([xmin,ymin]), np.array([xmax,ymax])]



    def closest_point(self, xx):
        """ set up Newton's method and initial guesses"""
        x, y = xx.T
        xs = self._xs
        ys = self._ys
        xp = self._xp
        yp = self._yp
        xpp = self._xpp
        ypp = self._ypp
        endpt1 = self._l
        endpt2 = self._r

        # distance squared from (x,y) to curve
        def d2(t, x, y):
            return (xs(t) - x)**2 + (ys(t) - y)**2

        # derivative of dist wrt t
        def g(t, x, y):
            return 2*(xs(t) - x)*xp(t) + 2*(ys(t) - y)*yp(t)

        # second derivative
        def gp(t, x, y):
            return 2*xp(t)*xp(t) + 2*(xs(t) - x)*xpp(t) + 2*yp(t)*yp(t) + 2*(ys(t) - y)*ypp(t)

        # Initial guess
        # use M equispaced samples in the parameter space
        M = 50
        s_guess = np.zeros_like(x)
        mindd_guess = np.zeros_like(x)
        ss = np.linspace(endpt1, endpt2, M, False)

        # for now, like cpParamCurveClosed_oldloop (matlab version)
        # ie process data all at once

        xlen = x.shape[0]

        # TODO: this is also written in a loop for now
        for it in xrange(xlen):
            dd = d2(ss, x[it], y[it])
            mindd_guess[it] = np.min(dd)
            i = np.argmin(dd)
            s_guess[it] = ss[i]
        
        cpx, cpy, dist, fail = self.newton(g, gp, x, y, xs, ys, s_guess, mindd_guess, endpt1, endpt2)

        cpxx = np.column_stack((cpx, cpy))
        return cpxx, dist, np.zeros(x.ndim), {}

    def newton(self, g, gp, x, y, xs, ys, s_guess, mindd_guess, endpt1, endpt2):
        """ Newton's method """
        tol = 1e-14
        maxn = 100;

        n = 1
        s = s_guess

        while True:
            n = n + 1

            # if the second-deriv of some components is close to zero, don't update s for those ones.  Otherwise do a newton step.
            snew = s
            I = 1*(abs(gp(s,x,y)) > tol)
            I = np.nonzero(I)[0]
            snew[I] = s[I] - g(s[I], x[I], y[I]) / gp(s[I], x[I], y[I])

            # Force parameter to be periodic
            snew2 = snew % (endpt2-endpt1) 
            #self.argPeriodic(snew-endpt1, endpt2-endpt1) + endpt1
            snew = snew2
            
            if n > maxn:
                fail = 1
                warnings.warn('max iterations')
                break
            elif max(abs(s-snew2)) < tol:
                fail = 0
                break
            
            # update
            s = snew;

        cpx = xs(s)
        cpy = ys(s)
        dist = np.sqrt((cpx - x)**2 + (cpy - y)**2)

        return cpx, cpy, dist, fail
    
    
    def argPeriodic(self, v, P=2*np.pi):
        """ Return the principal value of the argument """
        u = v % P
        return u


    def parametric_grid(self, rez=256):
        """
        Parameterized form (for plotting)
        """
        s = np.linspace(self._l, self._r, rez)
        X = self._xs(s)
        Y = self._ys(s)
        # TODO: return a list?
        return X, Y

