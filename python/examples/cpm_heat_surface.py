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
from petsc4py import PETSc

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
#    gv.view()
#    v = band.getCoordinates()
#    mlab.points3d(v[:,0],v[:,1],v[:,2], gv.getArray() ,mode = 'point')
#    band.g2l.scatter(gv,lv,PETSc.InsertMode.INSERT)
#    band.l2g.scatter(lv,wv,PETSc.InsertMode.INSERT)
#    mlab.figure()
#    mlab.points3d(v[:,0],v[:,1],v[:,2], wv.getArray() ,mode = 'point')
#    M.mult(wv,gv)
#    mlab.figure()
#    mlab.points3d(v[:,0],v[:,1],v[:,2], gv.getArray() ,mode = 'point')
    ts = PETSc.TS().create(comm=comm)
    ts.setProblemType(ts.ProblemType.LINEAR)
    ts.setType(ts.Type.EULER)  
#    mlab.show()
    ts.setTime(0.0)
    ts.setTimeStep(band.dx**2)
    ts.setMaxTime(1)
    ts.setMaxSteps(100)
    
    