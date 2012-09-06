'''
Created on Aug 10, 2012

@author: nullas
'''
import scipy as sp
try:
    from mayavi import mlab as pl
except ImportError:
    from enthought.mayavi import mlab as pl
from mpi4py import MPI
from cp.surfaces.MeshWrapper import MeshWrapper
from cp.petsc.band import Band
import sys
import petsc4py
from petsc4py import PETSc


petsc4py.init(sys.argv)


def test_initialu(cp):
    return sp.ones(cp.shape[0])#cp[:,0]
def initialu(cp):
    rlt = sp.ones(cp.shape[0])
    (ind,) = sp.where(2.3*cp[:,0]+cp[:,2] > 0)
    rlt[ind] = 0
    return rlt

def plot3d_sep(x,gv):
    return
    v2 = gv.getArray()
    pl.figure(bgcolor=(1,1,1),fgcolor=(0.5,0.5,0.5))
    pl.points3d(x[:,0],x[:,1],x[2::3],v2,mode='point')

def plot3d_total(x,gv = None):
    return
    v = PETSc.Vec().createWithArray(x)
    v = band.toZeroStatic(v)
    v = v.getArray()
    if gv is not None:
        v2 = band.toZeroStatic(gv)
        if comm.rank == 0:
            pl.figure(bgcolor=(1,1,1),fgcolor=(0.5,0.5,0.5))
            pl.points3d(v[0::3],v[1::3],v[2::3],v2.getArray(),mode='point')
    else:
        if comm.rank == 0:
            pl.figure(bgcolor=(1,1,1),fgcolor=(0.5,0.5,0.5))
            pl.points3d(v[0::3],v[1::3],v[2::3],mode='point')        
    
def plot3dformesh(x,cv,f):
    return
    cv = band.toZeroStatic(cv)
    if MPI.COMM_WORLD.rank == 0:
        v2 = cv.getArray()
        pl.figure(bgcolor=(1,1,1),fgcolor=(0.5,0.5,0.5))
        pl.triangular_mesh(x[:,0],x[:,1],x[:,2],f,scalars=v2)
        
def outputBin(gv):
    cv = band.toZeroStatic(gv)
    if comm.rank == 0:
        viewer = PETSc.Viewer().createBinary('gv.dat',mode='w',comm = MPI.COMM_SELF)
        cv.view(viewer)

if __name__ == '__main__':
    opt = {'M':100,'m':5,'d':3}
    try:
        surface = MeshWrapper('eight.ply')
        #surface = MeshWrapper('eight_refined.ply')
    except IOError:
        surface = MeshWrapper('cp/tests/data/eight.ply')
        #surface = MeshWrapper('cp/tests/data/eight_refined.ply')
    comm = MPI.COMM_WORLD
    band = Band(surface,comm,opt)
    la,lv,gv,wv = band.createGLVectors()
    v = band.getCoordinates()
    centers = band.BlockInd2CenterCarWithoutBand(band.gindBlockWBand.getArray())
    dt = 0.1*band.dx**2
    vv = sp.array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
    weights = sp.array([-6,1,1,1,1,1,1])*(dt/band.dx**2)
    L = band.createAnyMat(vv, weights, (7,3))
    PETSc.Sys.Print('Laplacian')

    M = band.createExtensionMatForLoop()
    band.initialu(initialu)
    PETSc.Sys.Print('Initial')
    plot3d_total(v,gv)
    nextt = 0.01
    PETSc.Sys.Print('Begin to solve.\n dt is {0}'.format(dt))
    for t in sp.arange(0,0.2,dt):
        L.multAdd(gv,gv,wv)
        M.mult(wv,gv)
        if t > nextt:
            nextt += 0.01
            PETSc.Sys.Print('time is {0}'.format(t))
            
    PETSc.Sys.Print('End to solve.')
    
    v = surface.v
    vsize = v.shape[0]
    vAssigned = vsize // comm.size + int(comm.rank < (vsize % comm.size))
    vstart = comm.exscan(vAssigned)
    PETSc.Sys.syncPrint(vAssigned)
    PETSc.Sys.syncFlush()
    if comm.rank == 0:
        vstart = 0
        
    mv = band.createExtensionMatForLoop(cp=v[vstart:vstart+vAssigned])
    PETSc.Sys.syncPrint('build extMat')
    PETSc.Sys.syncFlush()
    cv = mv.getVecLeft()
    mv.mult(gv,cv)
    outputBin(cv)
    PETSc.Sys.Print('maximal is {0}'.format(gv.max()[1]))    
del band,v   
    
    
    
    
