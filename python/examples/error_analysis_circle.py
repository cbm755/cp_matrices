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
from scipy.linalg import norm

petsc4py.init(sys.argv)


def test_initialu(cp):
    return sp.ones(cp.shape[0])#cp[:,0]
def initialu(cp):
    return cp[:,0]

def triplot(x,y,z,r=0.0002,title = 'band'):
#    z = c-c.min()
#    z /= z.max()
    return 0
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
    MBlocklist = [20,40,80,160]
    error = []
    dx = []
    l = sp.linspace(-sp.pi, sp.pi, 500)
    points = sp.column_stack((sp.cos(l),sp.sin(l)))
    exactu = sp.cos(l)
    
    comm = MPI.COMM_WORLD
    
    
    vsize = points.shape[0]
    vAssigned = vsize // comm.size + int(comm.rank < (vsize % comm.size))


    vstart = comm.exscan(vAssigned)
    if comm.rank == 0:
        vstart = 0
    for MBlock in MBlocklist:
        opt = {'M':MBlock,'m':5,'d':2}
        surface = Sphere(center=sp.array([0.0, 0.0]))
        
        band = Band(surface,comm,opt)
        la,lv,gv,wv = band.createGLVectors()
        v = band.getCoordinates() 
        dt = 0.1*band.dx**2
        vv = sp.array([[0,0],[1,0],[-1,0],[0,1],[0,-1]])
        weights = sp.array([-4,1,1,1,1])*(dt/band.dx**2)
        L = band.createAnyMat(vv, weights, (5,2))
        PETSc.Sys.Print('Laplacian')
    
        
        M = band.createExtensionMat()
        PETSc.Sys.Print('ExtensionMat built')
    
        band.initialu(initialu)
        PETSc.Sys.Print('Initial')
        nextt = 0.1
        PETSc.Sys.Print('Begin to solve')
        t = 0
        for t in sp.arange(0,1,dt):
            L.multAdd(gv,gv,wv)
            M.mult(wv,gv)
            if t > nextt:
                nextt += 0.1
                PETSc.Sys.Print('time is {0}'.format(t))
                

        mv = band.createExtensionMatForLoop(cp=points[vstart:vstart+vAssigned])
        cv = mv.getVecLeft()
        mv.mult(gv,cv)
        cv = band.toZeroStatic(cv)
        if comm.rank == 0:
            cv = cv.getArray()
            exu = sp.exp(-t)*exactu
            ee = norm(cv-exu, sp.inf)
            error.append(ee)
            dx.append( band.dx )
        
    
        
        PETSc.Sys.Print('==================================')   
        if comm.rank == 0: 
            print('maximal is {0}'.format(ee))
        PETSc.Sys.Print('==================================')   
        del band,v,mv,cv 
    if comm.rank == 0:
        import pickle
        ferror = open('error.pickle','w')
        fdx = open('dx.pickle','w')
        pickle.dump(error, ferror)
        pickle.dump(dx,fdx)
    
    
    