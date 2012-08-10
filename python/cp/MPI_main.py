'''
Created on Jul 29, 2012

@author: nullas
'''
from __future__ import division
from surfaces import Sphere as Surface
from petsc.band import Band
from mpi4py import MPI


try:
    import petsc4py
    import sys
    petsc4py.init(sys.argv)
except Exception as exp:
    print exp
    


if __name__ == '__main__':
    surface = Surface()
    bnd = Band(surface,MPI.COMM_WORLD)
    bnd.createGLVectors()
    M = bnd.createExtensionMat()
    