"""
A class to make a new closest point representation from the union of a
list of other closest point representations.

Algorithm: take the minimum distance and return that cp

TODO: dealing with boundaries properly might be difficult
"""
from Surface import Surface
#import numpy
from numpy import array as a

class Multi(Surface):
    def __init__(self, list):
        self._list = list

    def closestPointToCartesian(self, x):
        print 'TODO: in multi, WIP'
        mindist = inf
        for l in list:
            cp,dist,bdy,other = l.closestPointToCartesian(x)
            if dist < mindist:
                mincp = cp
                mindist = dist
                minbdy = bdy
                min_other = other
                whichPiece = 66   # todo, index in list
        min_other['whichPiece'] = whichPiece

        return mincp, mindist, minbdy, min_other
