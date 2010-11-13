"""
CPBar function.  Wraps a object and does the CPBar projection to it if
the closest point hit a bdy.
"""
from ClosestPoint import ClosestPoint

import numpy
from numpy import array as a

class CPBar(ClosestPoint):
    def __init__(self, parent):
        self._obj = parent

    def closestPointToCartesian(self, x):
        cpx,dist,bdy = self._obj.closestPointToCartesian(x)
        if bdy==1 or bdy==2:
            y = x + 2*(cpx - x)
            cpx2,dist2,bdy2 = cpf(y)
            if (bdy2 != 0):
                #raise NameError('cpbar hit bdy!  What to do?')
                print 'cpbar hit bdy!  dist=',dist
                print x
                print dist,cpx,bdy
                print dist2,cpx2,bdy2
        else:
            cpx2 = cpx
        return (cpx2, dist, bdy)
