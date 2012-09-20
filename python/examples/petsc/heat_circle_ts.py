'''
Created on Aug 10, 2012

@author: nullas
'''
import numpy as np
import timeit
from matplotlib import pylab as pl
from mpi4py import MPI
#from cp.surfaces.MeshWrapper import MeshWrapper
from cp.surfaces.Sphere import Sphere
from cp.surfaces.coordinate_transform import cart2pol
from cp.petsc.band_with_doc import Band
import sys
import petsc4py
from petsc4py import PETSc
import matplotlib.tri as tri
from scipy.linalg import norm

petsc4py.init(sys.argv)


def uexactfn(t,cp):
    k1 = 2
    k2 = 4
    k3 = 6
    th, r = cart2pol(cp[:,0],cp[:,1])
    return np.exp(-k1**2*t)*np.cos(k1*th) +  np.exp(-k2**2*t)*np.sin(k2*th)\
           + np.exp(-k3**2*t)*np.cos(k3*th)

def uinitialfn(cp):
    return uexactfn(0,cp)

def monitor(ts,i,t,x):
    PETSc.Sys.Print(
        'At {0} max-norm is {1}'.format(t, x.norm(PETSc.NormType.INFINITY)))


if __name__ == '__main__':
    MBlocklist = [20,40,80,160]
    error = []
    dx = []
    l = np.linspace(-np.pi, np.pi, 500)
    points = np.column_stack((np.cos(l),np.sin(l)))
    exactu = np.cos(l)
    
    comm = MPI.COMM_WORLD
   
    # these are for plotting 
    vsize = points.shape[0]
    vAssigned = vsize // comm.size + int(comm.rank < (vsize % comm.size))
    vstart = comm.exscan(vAssigned)
    if comm.rank == 0:
        vstart = 0

    for MBlock in MBlocklist:
        opt = {'M':MBlock,'m':5,'d':2}
        surface = Sphere(center=np.array([0.0, 0.0]))
        
        band = Band(surface,comm,opt)
        _,_,v,v2 = band.createGLVectors()
        grid = band.getCoordinates() 
        band.computeCP()

        #vv = np.array([[0,0],[1,0],[-1,0],[0,1],[0,-1]])
        #weights = np.array([-4,1,1,1,1])*(dt/band.dx**2)
        #L = band.createAnyMat(vv, weights, (5,2))
        L = band.createLaplacianMat()
        PETSc.Sys.Print('Laplacian built')
    
        
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
        Tf = 1
        dt = 0.2*band.dx**2
        numtimesteps = int( Tf // dt + 1 )

        ts = PETSc.TS().create(comm=comm)
        ts.setProblemType(ts.ProblemType.LINEAR)
        ts.setType('euler')
#        ts.setType('beuler')
        #ts.setType('cn')  # Crank-Nicolson, very slow
#        ts.setType('ssp')
        # Lots of printout. How to turn this off?
#        ts.setType('rk')
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
        ksp.setType(ksp.Type.GMRES)        # GMRES
        pc = ksp.getPC()  
#        pc.setType('none')

        PETSc.Sys.Print('Begin to solve')
        start_time = timeit.default_timer()
        t = ts.solve(v)

        print "Times, rank=", comm.rank, "time=", timeit.default_timer() - start_time

        #mv = band.createExtensionMatForLoop(cp=points[vstart:vstart+vAssigned])
        Eplot = band.createExtensionMat(cp=points[vstart:vstart+vAssigned])
        cv = Eplot.getVecLeft()
        Eplot.mult(v,cv)
        cv = band.toZeroStatic(cv)
        if comm.rank == 0:
            cv = cv.getArray()
            exu = uexactfn(t,points)
            ee = norm(cv-exu, np.inf) / norm(exu, np.inf)
            error.append(ee)
            dx.append( band.dx )
        
        PETSc.Sys.Print('==================================')   
        if comm.rank == 0: 
            print('final time is {0}'.format(t))
            print('maximal is {0}'.format(ee))
        PETSc.Sys.Print('==================================')   
        del band,v,Eplot,cv 
    if comm.rank == 0:
        import pickle
        ferror = open('error.pickle','w')
        fdx = open('dx.pickle','w')
        pickle.dump(error, ferror)
        pickle.dump(dx,fdx)
    
    
    
