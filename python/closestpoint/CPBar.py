"""
CPBar function.  Wraps a object and does the CPBar projection to it if
the closest point hit a bdy.
"""
from ClosestPoint import ShapeWithBdy

import numpy
from numpy import array as a

class CPBar(ShapeWithBdy):
    def __init__(self, parent):
        self._obj = parent
        self._dim = self._obj._dim
        self._hasParam = self._obj._hasParam
        #self._bb = self._obj

    def ParamGrid(self,rez=None):
        if rez==None:
            return self._obj.ParamGrid()
        else:
            return self._obj.ParamGrid(rez=rez)
    def getBB(self):
        return self._obj.getBB()

    def closestPointToCartesian(self, x):
        cpf = self._obj.closestPointToCartesian
        cpx,dist,bdy = cpf(x)
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
    cp = closestPointToCartesian
