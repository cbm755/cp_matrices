"""
Tools for manipulating closest point representations.

Currently does translation (shift).

TODO: implement rotation, use a matrix internally.

TODO: should deal with bdy, others in a reasonable way...
"""
from Surface import Surface
from numpy import array as a

class Translate(Surface):
    def __init__(self, cpparent, shift):
        self._shift = shift
        self._cpparent = cpparent

    def closestPointToCartesian(self, x):
        #cp,dist = super(x - self._shift)
        # TODO: need to think about how to deal with bdy in this case,
        # for now just pass it onwards.
        cp,dist,bdy,others = self._cpparent.closestPointToCartesian(x - self._shift)
        cp = cp + self._shift
        return cp, dist, bdy, others
    cp = closestPointToCartesian

    # TODO: vararg here
    def ParamGrid(self):
        """
        TODO: currently hardcoded for 3D, and single component too
        """
        sh = self._shift
        # TODO: should pass varg here
        x,y,z = self._cpparent.ParamGrid()
        x,y,z = x,y,z + sh[0],sh[1],sh[2]
        return x,y,z
