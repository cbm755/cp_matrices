"""Solves the heat equation on a true sphere."""

import numpy as np
from matplotlib import pylab as pl
from mpi4py import MPI
from cp.surfaces import Sphere
from cp.petsc.band_with_doc import Band
import sys
import petsc4py
from petsc4py import PETSc
from scipy.linalg import norm

try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab

from cp.surfaces.coordinate_transform import cart2sph

def initialu(cp):
    th,phi,r = cart2sph(cp[:,0],cp[:,1],cp[:,2])
    return np.cos(phi + np.pi / 2)

if __name__ == '__main__':
    MBlocklist = [10,20,40,80]
    Tf = 0.5
    error = []
    dx = []
    
    comm = MPI.COMM_WORLD
    
    surface = Sphere(center=np.array([0.0, 0.0, 0.0]))
    # exact solution  
    rez = 30
    xx, yy, zz = surface.parametric_grid(rez)
    th, phi, r = cart2sph(xx,yy,zz)
    exactu = np.cos(phi + np.pi / 2).ravel()
    
    points = np.array([xx.ravel(),yy.ravel(),zz.ravel()]).T
     
    vsize = xx.ravel().shape[0]
    vAssigned = vsize // comm.size + int(comm.rank < (vsize % comm.size))

    vstart = comm.exscan(vAssigned)
    if comm.rank == 0:
        vstart = 0
    for MBlock in MBlocklist:
        opt = {'M':MBlock,'m':5,'d':3}

        
        band = Band(surface,comm,opt)
        la,lv,gv,wv = band.createGLVectors()
        grid = band.getCoordinates() 
        band.computeCP()


        dt = 0.2*band.dx**2
        #vv = sp.array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
        #weights = sp.array([-6,1,1,1,1,1,1])*(dt/band.dx**2)
        #L = band.createAnyMat(vv, weights, (7,3))
        L = band.createLaplacianMat()
        L *= dt

        PETSc.Sys.Print('Laplacian matrix built')
        
        M = band.createExtensionMat()
        PETSc.Sys.Print('ExtensionMat built')
    
        # initial conditions, in serial code we could do:
        # th, phi, r = cart2sph(grid[:, 0], grid[:, 1], grid[:, 2])
        # u = np.cos(phi + np.pi / 2)
        # but in parallel code, we might need to do:
        gv.setArray(initialu(band.cp)) 
        PETSc.Sys.Print('Initial vector has been set up')
        nextt = 0.1
        PETSc.Sys.Print('Begin to solve')
        t = 0
        for t in np.arange(0,Tf,dt):
            L.multAdd(gv,gv,wv)
            M.mult(wv,gv)
            if t > nextt:
                nextt += 0.1
                PETSc.Sys.Print('time is {0}'.format(t))
                

        #mv = band.createExtensionMatForLoop(cp=points[vstart:vstart+vAssigned])
        mv = band.createExtensionMat(cp=points[vstart:vstart+vAssigned])
        cv = mv.getVecLeft()
        mv.mult(gv,cv)
        cv = band.toZeroStatic(cv)
        if comm.rank == 0:
            cv = cv.getArray()
            exu = np.exp(-2*t)*exactu
            ee = norm(cv-exu, np.inf)
            error.append(ee)
            dx.append( band.dx )
        
    
        
        PETSc.Sys.Print('==================================')   
        if comm.rank == 0: 
            print('maximal is {0}'.format(ee))
        PETSc.Sys.Print('==================================')   
        del band,grid,mv,cv 
    if comm.rank == 0:
        import pickle
        ferror = open('error.pickle','w')
        fdx = open('dx.pickle','w')
        pickle.dump(error, ferror)
        pickle.dump(dx,fdx)
    
    
    



