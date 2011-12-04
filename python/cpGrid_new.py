# make print a function as in python 3
#from __future__ import print_function

import numpy
#from numpy import array as numpy.array

from cpOps import findGridInterpBasePt, buildInterpWeights

class cpGrid:
    """
    store a cp-matrix-style grid
    """
    def __init__(self):
        0

    def __init__(self, x1d, y1d, dx):
        self.x1d = x1d
        self.y1d = y1d
        self.dx = dx
        self.nx = x1d.shape[0]
        self.ny = y1d.shape[0]
        #self.level = level
        #self.gridpt = lsc
        #self.dx = dx
        #self.basicPt = basicPt
        #self.children = []

    def __str__(self):
        return 'TODO'
        #return 'gridpt ' + str(self.gridpt) + \
        #    ' (level ' + str(self.level) + ')'

    def sub2ind(self, sub):
        return s2i( (self.nx,self.ny), sub)
    
    def ind2sub(self, ind):
        return i2s( (self.nx,self.ny), ind)

def i2s(sz,ind):
    """works like ind2sub but with v = s2i(sz,ind) instead of
    [v1,v2,v3,...] = ind2sub(sz,ind)
    works along dimension 2 (i.e. rows specify the index)
    CAUTION no error checking (if index out of range)"""

    ind = numpy.array(ind) + 1
    sz = numpy.array(sz)
    sub=numpy.zeros([ind.size,sz.size]).astype(numpy.int);

    for i in range(len(sz)-1):

        sub[:,i] = numpy.mod(ind,sz[i]);
        sub[sub[:,i]==0,i] = sz[i];
        ind = (ind - sub[:,i])/sz[i] + 1;

    sub[:,i+1] = ind;
    sub = sub - 1

    return sub



def s2i(sz,v):
    """ works like sub2ind but with s2i(sz,[v1,v2,v3,...]) instead of
    sub2ind(sz,v1,v2,v3,...)
    works along dimension 2 (i.e. rows specify the indices)
    CAUTION no error checking (if index out of range)"""

    index = numpy.array(v[:,0]).astype(numpy.int)
    sc=1
    for i in range(1,len(sz)):
        sc = sc*(sz[i-1]);
        index = index + sc*(v[:,i]);

    return index


if __name__ == "__main__":
    print "running as a script, import it instead"
