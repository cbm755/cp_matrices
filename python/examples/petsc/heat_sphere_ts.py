"""Solves the heat equation on a true sphere."""

import numpy as np
import timeit
#from matplotlib import pylab as pl
from mpi4py import MPI
from cp.surfaces import Sphere
from cp.petsc.band_with_doc import Band
import sys
import petsc4py
from petsc4py import PETSc
from scipy.linalg import norm

#try:
#    from mayavi import mlab
#except ImportError:
#    from enthought.mayavi import mlab

from cp.surfaces.coordinate_transform import cart2sph

def uexactfn(t,cp):
    th,phi,r = cart2sph(cp[:,0],cp[:,1],cp[:,2])
    return np.exp(-2*t)*np.cos(phi + np.pi / 2)

if __name__ == '__main__':
    MBlocklist = [10,20,40,80]
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
        _,_,v,v2 = band.createGLVectors()
        grid = band.getCoordinates() 
        band.computeCP()

        #vv = sp.array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
        #weights = sp.array([-6,1,1,1,1,1,1])*(dt/band.dx**2)
        #L = band.createAnyMat(vv, weights, (7,3))
        L = band.createLaplacianMat()
        PETSc.Sys.Print('Laplacian matrix built')

        E = band.createExtensionMat()
        E1 = band.createExtensionMat(p=1)
        PETSc.Sys.Print('ExtensionMat built')
    
        v.setArray(uexactfn(0,band.cp))
        PETSc.Sys.Print('Initial value has been set')
        
        la = 2*band.Dim / band.dx**2
        # TODO: easier way to set up the identity matrix?
        I = PETSc.Mat().create(comm=comm)
        I.setSizes((v.sizes,v.sizes))
        I.setFromOptions()
        I.setPreallocationNNZ(1)
        v2.set(1)
        I.setDiagonal(v2)
        I.assemble()
            
        # TODO: clearer way to compute M = EL-lambda*(I-E)?
        M = E1.matMult(L)
        IminusE = I.copy()
        IminusE.axpy(-1,E)
        M.axpy(-la,IminusE)


        t = 0
        Tf = 0.2
        dt = 0.2*band.dx**2
        numtimesteps = int( Tf // dt + 1 )

        ts = PETSc.TS().create(comm=comm)
        ts.setProblemType(ts.ProblemType.LINEAR)
        ts.setType('euler') # Forward Euler

#        dt = 0.2*band.dx
#        ts.setType('beuler') # Backward Euler, diverge, don't know why

#        ts.setType('ssp')  # slower than euler, more accurate
#        ts.setType('cn')   # Crank-Nicolson, very slow
        # Lots of printout. How to turn this off?
#        ts.setType('rk')
#        ts.setType('arkimex') # diverge

	ts.setSolution(v)
        ts.setFromOptions()
        ts.setRHSFunction(PETSc.TS.computeRHSFunctionLinear,v2)
        ts.setRHSJacobian(PETSc.TS.computeRHSJacobianConstant,M,M)
        #ts.setMonitor(monitor)
        ts.setTime(0.0)
        ts.setInitialTimeStep(0.0,dt)
        ts.setMaxSteps(numtimesteps)
        ts.setMaxSNESFailures(-1)       # allow an unlimited number of failures (step will be rejected and retried)

        snes = ts.getSNES()             # Nonlinear solver
        ksp = snes.getKSP()             # Linear solver
        ksp.setType(ksp.Type.GMRES)     # GMRES
#        ksp.setType('bcgs')             # Bicgstab
        ksp.setTolerances(atol=1e-10,rtol=1e-10)
        pc = ksp.getPC()  
        pc.setType('none')
#        pc.setType('gamg')

        PETSc.Sys.Print('Begin to solve')
        start_time = timeit.default_timer()
        t = ts.solve(v)

        print "Times, rank=", comm.rank, "time=", timeit.default_timer() - start_time

        #mv = band.createExtensionMatForLoop(cp=points[vstart:vstart+vAssigned])
        mv = band.createExtensionMat(cp=points[vstart:vstart+vAssigned])
        cv = mv.getVecLeft()
        mv.mult(v,cv)
        cv = band.toZeroStatic(cv)
        if comm.rank == 0:
            cv = cv.getArray()
            exu = np.exp(-2*t)*exactu
            ee = norm(cv-exu, np.inf) / norm(exu, np.inf)
            error.append(ee)
            dx.append( band.dx )
        
        
        PETSc.Sys.Print('==================================')   
        if comm.rank == 0: 
            print('final time is {0}'.format(t))
            print('maximal is {0}'.format(ee))
        PETSc.Sys.Print('==================================')   
        del band,grid,mv,cv 
    if comm.rank == 0:
        import pickle
        ferror = open('error.pickle','w')
        fdx = open('dx.pickle','w')
        pickle.dump(error, ferror)
        pickle.dump(dx,fdx)
    
    
    



