'''
Created on Aug 10, 2012

@author: nullas
'''
import scipy as sp
from mayavi import mlab as pl
from mpi4py import MPI
from cp.surfaces.MeshWrapper import MeshWrapper
from cp.surfaces.Sphere import Sphere
from cp.petsc.band import Band
import sys
import petsc4py
from petsc4py import PETSc


petsc4py.init(sys.argv)


def test_initialu(cp):
    return sp.ones(cp.shape[0])#cp[:,0]
def initialu(cp):
    return cp[:,0]

def plot3d(x,gv):
    v1 = PETSc.Vec().createWithArray(x)
    v = band.toZeroStatic(v1)
    v2 = band.toZeroStatic(gv)
    pl.figure()
    pl.points3d(v[0::3],v[1::3],v[2::3],v2,mode='point')
    
def plot3dformesh(x,gv):
    v2 = band.toZeroStatic(gv)
    if MPI.COMM_WORLD.rank == 1:
        pl.figure()
        pl.points3d(x[:,0],x[:,1],v[:,2],v2,mode='point')

if __name__ == '__main__':
    opt = {'M':80,'m':5,'d':3}
    surface = MeshWrapper('knot1.ply')
    pl.triangular_mesh(surface.v[:,0],surface.v[:,1],surface.v[:,2],surface.f,opacity = 1)
    comm = MPI.COMM_WORLD
    band = Band(surface,comm,opt)
    la,lv,gv,wv = band.createGLVectors()
    v = band.getCoordinates() 
    pl.points3d(v[:,0],v[:,1],v[:,2],mode='point')
    dt = 0.3*band.dx**2
    vv = sp.array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
    weights = sp.array([-6,1,1,1,1,1,1])*(dt/band.dx**2)
    L = band.createAnyMat(vv, weights, (7,3))
    PETSc.Sys.Print('Laplacian')


    
    M = band.createExtensionMatForLoop()

    


#    LM = M.copy()
#    M.matMult(L,LM)
#    ts = PETSc.TS().create(comm=comm)
#    ts.setProblemType(ts.ProblemType.LINEAR)
#    ts.setType(ts.Type.EULER)  
#    ts.setTime(0.0)
#    ts.setTimeStep(band.dx**2)
#    ts.setMaxTime(1)
#    ts.setMaxSteps(1000)
#    ts.setSolution(gv)
#    ts.setFromOptions()
#    ts.setRHSFunction(None,wv)
#    ts.setRHSJacobian(None,LM,LM)
#    ts.solve(gv)
    band.initialu(initialu)
    PETSc.Sys.Print('Initial')
    plot3d(v,gv)
    nextt = 0.1
    PETSc.Sys.Print('Begin to solve')
    for t in sp.arange(0,0,dt):
        L.multAdd(gv,gv,wv)
        M.mult(wv,gv)
        if t > nextt:
            nextt += 0.1
            PETSc.Sys.Print('time is {0}'.format(t))



    plot3dformesh(v,gv)
    pl.show()
    PETSc.Sys.Print('maximal is {0}'.format(gv.max()[1]))    
del band,v   
    
    
    
    