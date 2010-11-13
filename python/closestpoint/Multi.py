"""
TODO
"""
from ClosestPoint import ClosestPoint

#import numpy
from numpy import array as a

class Multi(ClosestPoint):
    def __init__(self, list):
        self._list = list

    def closestPointToCartesian(self, x):
        print 'hello in multi, WIP'
        mindist = inf
        for l in list:
            cp,dist = l.closestPointToCartesian(x)
            if dist < mindist:
                mincp = cp
                mindist = dist
                whichPiece = 66   # todo
        return mincp,mindist,whichPiece
