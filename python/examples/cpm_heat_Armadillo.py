'''
Created on Aug 10, 2012

@author: nullas
'''
import scipy as sp
from mayavi import mlab
from mpi4py import MPI
from cp.surfaces.MeshWrapper import MeshWrapper
from cp.surfaces.Sphere import Sphere
from cp.petsc.band import Band
import sys
import petsc4py
from petsc4py import PETSc

petsc4py.init(sys.argv)
def initialu(cp):
    return cp[:,0]**2+10*cp[:,1]



if __name__ == '__main__':
    opt = {'M':40,'m':5}
    surface = Sphere()
#    mlab.triangular_mesh(surface.v[:,0],surface.v[:,1],surface.v[:,2],surface.f,opacity = 0.2)
    comm = MPI.COMM_WORLD
    band = Band(surface,comm,opt)
    la,lv,gv,wv = band.createGLVectors()
    M = band.createExtensionMat()
    PETSc.Sys.Print('ExtensionMat built')
    band.initialu(initialu)
    PETSc.Sys.Print('Initial')
    dt = 0.3*band.dx**2
    v = sp.array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
    weights = sp.array([-6+dt/band.dx**2,1,1,1,1,1,1])*(dt/band.dx**2)
    L = band.createAnyMat(v, weights, (7,3))
    PETSc.Sys.Print('Laplacian')
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

    
    nextt = 0.01
    PETSc.Sys.Print('Begin to solve')
    for t in sp.arange(0,0.1,dt):
        L.mult(gv,wv)
        M.mult(wv,gv)
        if t > nextt:
            nextt += 0.01
            PETSc.Sys.Print('time is {0}'.format(t))
#    L.mult(gv,wv)
#    M.mult(wv,gv)    
    v = band.getCoordinates() 
    mlab.points3d(v[:,0],v[:,1],v[:,2],wv.getArray(),mode = 'point')
    mlab.show()
#    wv.view()
    
del band,v   
    
    