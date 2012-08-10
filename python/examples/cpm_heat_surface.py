'''
Created on Aug 10, 2012

@author: nullas
'''
import scipy as sp
from mayavi import mlab
from mpi4py import MPI
from cp.surfaces.MeshWrapper import MeshWrapper
from cp.petsc.band import Band
from cp.surfaces.Sphere import Sphere

def initialu(cp):
    return sp.cos(cp[:,2])+cp[:,0]


if __name__ == '__main__':
    opt = {'M':40,'m':5}
    surface = Sphere()
    comm = MPI.COMM_WORLD
    band = Band(surface,comm,opt)
    la,lv,gv,wv = band.createGLVectors()
    band.getCP()
    M = band.createExtensionMat()
    band.initialu(initialu)
    gv.view()
    v = band.getCoordinates()
    mlab.points3d(v[:,0],v[:,1],v[:,2], gv.getArray() ,mode = 'point')
    mlab.show()
    band.g2l.scatter(gv,lv,PETSc.InsertMode.INSERT)
    
    