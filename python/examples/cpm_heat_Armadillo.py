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

def plot3d(x,c,title = 'band'):
#    z = c-c.min()
#    z /= z.max()
    pl.figure()
    pl.points3d(x[:,0],x[:,1],x[:,2],c,mode='point')

if __name__ == '__main__':
    opt = {'M':50,'m':10,'d':3}
    surface = MeshWrapper()
    pl.triangular_mesh(surface.v[:,0],surface.v[:,1],surface.v[:,2],surface.f,opacity = 0.2)
    comm = MPI.COMM_WORLD
    band = Band(surface,comm,opt)
    la,lv,gv,wv = band.createGLVectors()
    v = band.getCoordinates() 
    dt = 0.3*band.dx**2
    vv = sp.array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
    weights = sp.array([-6,1,1,1,1,1,1])*(dt/band.dx**2)
    L = band.createAnyMat(vv, weights, (7,3))
    PETSc.Sys.Print('Laplacian')
    band.test_initialu(test_initialu)
    L.mult(gv,wv)
    c = wv.getArray()
#    wv.view()
    plot3d(v,c,title='Lu')
#    wv.view()
    
    M = band.createExtensionMat()
    PETSc.Sys.Print('ExtensionMat built')
#    band.initialu(initialu)
    
    M.mult(gv,wv)
    c = wv.getArray()
    plot3d(v,c,title='Mu')
#    (ind,) = sp.where(c>1.000000000001)
#    print 'nonzero terms>>1'
#    for i in ind:
#        b = M[i,:]
#        print i
#        print b[b.nonzero()[0]]
#        print M[i,:].sum()
    
#    (ind,) = sp.where(c<0.999999999999)
#    print 'nonzero terms<<1'
#    for i in ind:
#        b = M[i,:]
#        print i
#        print b[b.nonzero()[0]]
#        print M[i,:].sum()

    c = gv.getArray()

#    c -= c.min()
#    c /= c.max()
#    gv.view()
#    pl.scatter(v[:,0],v[:,1],c=c)
#    pl.axis('equal')
    plot3d(v,c,title='initial u')

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
    c = gv.getArray()
    plot3d(v,c,title='initial u')
    nextt = 0.1
    PETSc.Sys.Print('Begin to solve')
    for t in sp.arange(0,0.01,dt):
        L.multAdd(gv,gv,wv)
        M.mult(wv,gv)
        if t > nextt:
            nextt += 0.1
            PETSc.Sys.Print('time is {0}'.format(t))
#    L.mult(gv,wv)
#    M.mult(wv,gv)    

#    mlab.points3d(v[:,0],v[:,1],wv.getArray(),mode = 'point')
#    mlab.show()
#    wv.view()
#    pl.figure()
    c = gv.getArray()
#    gv.view()
#    c /= c.max()
#    pl.scatter(v[:,0],v[:,1],c=c)
#    pl.axis('equal')

    plot3d(v,c,title = 'Result')
    pl.show()
    PETSc.Sys.Print('maximal is {0}'.format(gv.max()[1]))    
del band,v   
    
    
    
    