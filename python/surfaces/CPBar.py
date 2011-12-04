"""
CPBar function.  Wraps a object and does the CPBar projection to it if
the closest point hit a bdy.
"""
from Surface import ShapeWithBdy

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
        """ note: returns distance to the original cp """
        #cpf = self._obj.closestPointToCartesian
        cpf = self._obj.cp
        cpx,dist,bdy,others = cpf(x)
        if bdy==1 or bdy==2:
            y = x + 2*(cpx - x)
            cpx2,dist2,bdy2,others2 = cpf(y)
            others['origcp'] = cpx
            others['others2'] = others2
            if (bdy2 != 0):
                #raise NameError('cpbar hit bdy!  What to do?')
                print 'cpbar hit bdy!  dist=',dist
                print x
                print dist,cpx,bdy
                print dist2,cpx2,bdy2
        else:
            cpx2 = cpx

        return cpx2, dist, bdy, others
    cp = closestPointToCartesian
