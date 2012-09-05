'''
Created on Aug 10, 2012

@author: nullas
'''
import scipy as sp
from matplotlib import pylab as pl
from mpi4py import MPI
#from cp.surfaces.MeshWrapper import MeshWrapper
from cp.surfaces.Sphere import Sphere
from cp.petsc.band import Band
import sys
import petsc4py
from petsc4py import PETSc
import matplotlib.tri as tri

petsc4py.init(sys.argv)
def test_initialu(cp):
    return sp.ones(cp.shape[0])#cp[:,0]
def initialu(cp):
    return cp[:,0]

def triplot(x,y,z,r=0.001,title = 'band'):
#    z = c-c.min()
#    z /= z.max()
    triang = tri.Triangulation(x,y)
    xmid = x[triang.triangles].var(axis=1)
    ymid = y[triang.triangles].var(axis=1)
    mask = sp.where(xmid*xmid + ymid*ymid > r*r, 1, 0)
    triang.set_mask(mask)
    pl.figure()
    pl.gca().set_aspect('equal')
    pl.tricontourf(triang, z)
    pl.colorbar()
    V = sp.arange(-10,10,dtype=sp.double)/10*z.max()
    pl.tricontour(triang, z,V)#, colors='k')
    pl.title(title)

if __name__ == '__main__':
    opt = {'M':40,'m':4,'d':2}
    surface = Sphere(center=sp.array([0.0, 0.0]))
    comm = MPI.COMM_WORLD
    band = Band(surface,comm,opt)
    la,lv,gv,wv = band.createGLVectors()
    v = band.getCoordinates() 
    dt = 0.1*band.dx**2
    vv = sp.array([[0,0],[1,0],[-1,0],[0,1],[0,-1]])
    weights = sp.array([-4,1,1,1,1])*(dt/band.dx**2)
    L = band.createAnyMat(vv, weights, (5,2))
    PETSc.Sys.Print('Laplacian')
    band.test_initialu(test_initialu)
    L.mult(gv,wv)
    c = wv.getArray()
    triplot(v[:,0],v[:,1],c,title='Lu')

    
    M = band.createExtensionMat()
    PETSc.Sys.Print('ExtensionMat built')

    
    M.mult(gv,wv)
    c = wv.getArray()
    triplot(v[:,0],v[:,1],c,title='Mu')


    c = gv.getArray()

    triplot(v[:,0],v[:,1],c,title='initial u')

    band.initialu(initialu)
    PETSc.Sys.Print('Initial')
    nextt = 0.1
    PETSc.Sys.Print('Begin to solve')
    for t in sp.arange(0,1,dt):
        L.multAdd(gv,gv,wv)
        M.mult(wv,gv)
        if t > nextt:
            nextt += 0.1
            PETSc.Sys.Print('time is {0}'.format(t))

    c = gv.getArray()


    triplot(v[:,0],v[:,1],c,title = 'Result')
    pl.show()
    PETSc.Sys.Print('maximal is {0}'.format(gv.max()[1]))    
del band,v   
    
    