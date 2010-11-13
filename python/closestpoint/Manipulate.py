"""
Tools for manipulating closest point representations
"""
#import numpy
from numpy import array as a


class Shift(ClosestPoint):
    def __init__(self, cpparent, shift):
        self._shift = shift
        self._cpparent = cpparent

    def closestPointToCartesian(self, x):
        #cp,dist = super(x - self._shift)
        cp,dist = self._cpparent.closestPointToCartesian(x - self._shift)
        cp = cp + self._shift
        return cp,dist

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
