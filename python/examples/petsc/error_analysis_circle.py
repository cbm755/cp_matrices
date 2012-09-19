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
    k = 2
    th, r = cart2pol(cp[:,0],cp[:,1])
    return np.exp(-k**2*t)*np.cos(k*th) 

def uinitialfn(cp):
    return uexactfn(0,cp)

def triplot(x,y,z,r=0.0002,title = 'band'):
#    z = c-c.min()
#    z /= z.max()
    return 0
    triang = tri.Triangulation(x,y)
    xmid = x[triang.triangles].var(axis=1)
    ymid = y[triang.triangles].var(axis=1)
    mask = np.where(xmid*xmid + ymid*ymid > r*r, 1, 0)
    triang.set_mask(mask)
    pl.figure()
    pl.gca().set_aspect('equal')
    pl.tricontourf(triang, z)
    pl.colorbar()
    V = np.arange(-10,10,dtype=np.double)/10*z.max()
    pl.tricontour(triang, z,V)#, colors='k')
    pl.title(title)

if __name__ == '__main__':
    MBlocklist = [20,40,80,160]
    error = []
    dx = []
    l = np.linspace(-np.pi, np.pi, 1000)
    points = np.column_stack((np.cos(l),np.sin(l)))
    
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
        
        # TODO: pass 'cpm' as a paramter when calling this example would be nicer;
        # currently we need to modify the code manually
        cpm = 3

        if cpm >= 1:  # Method of lines, von Glehn--Maerz--Macdonald
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
        dt = 0.1*band.dx**2
        numtimesteps = int( Tf // dt + 1 )
        nextt = 0.1

        PETSc.Sys.Print('Begin to solve')
        start_time = timeit.default_timer()

        # various CPM algorithms, if cpm>=1: use various MOLs, vGMM
        if cpm == 0:  # explicit euler, ruuth--merriman
            L *= dt
            for kt in xrange(numtimesteps):
                L.multAdd(v,v,v2)  # v2 = v + (dt*L)*v
                E.mult(v2,v)       # v = E*v2
                t = kt*dt
                if t > nextt:
                    nextt += 0.1
                    PETSc.Sys.Print('time is {0}'.format(t))

        elif cpm == 1:  # explicit Euler, von Glehn--Maerz--Macdonald
            M.scale(dt)
            for kt in xrange(numtimesteps):
                M.multAdd(v, v, v2)    # v2 = v + (dt*M)*v
                # todo: how to just assign data in v2 to v?
                v2.swap(v)
                t = kt * dt
                if t > nextt:
                    nextt += 0.1
                    PETSc.Sys.Print('time is {0}'.format(t))

        elif cpm == 2:   # implicit Euler
            dt = 0.1*band.dx
            numtimesteps = int( Tf // dt + 1 )
            A = I.copy()
            A.axpy(-dt,M)
            A.setUp()
            ksp = PETSc.KSP().create()
            ksp.setOperators(A)
            ksp.setType(ksp.Type.GMRES)
            pc = ksp.getPC()
            # slow
            #pc.setType(pc.Type.GAMG)
            # default preconditioner, fastest of all , and seems to be enough accurate...
            pc.setType('none')
            for kt in xrange(numtimesteps):
                ksp.solve(v,v2)
                v = v2
                t = kt*dt
                if t > nextt:
                    nextt += 0.1
                    PETSc.Sys.Print('time is {0}'.format(t))

        elif cpm == 3:   # Crank-Nicolson
            # Currently only first order, do not know why..
            dt = 0.1*band.dx
            numtimesteps = int( Tf // dt + 1 )
            beta = 0.5
            #A = I.copy()
            #A.axpy(-beta*dt,M)
            A = I - beta*dt*M
            A.setUp()
            ksp = PETSc.KSP().create()
            ksp.setOperators(A)
            #ksp.setType(ksp.Type.GMRES)
            ksp.setType(ksp.Type.PREONLY)   # Just use the preconditioner as a direct method
            pc = ksp.getPC()
            pc.setType(pc.Type.LU)          # This only works for a single processor 
            # slow
            #pc.setType(pc.Type.GAMG)
            # default preconditioner, fastest of all , and seems to be enough accurate...
            #pc.setType('none')
            M.scale((1.-beta)*dt)
            for kt in xrange(numtimesteps):
                M.multAdd(v, v, v2)    # v2 = v + (0.5*dt*M)*v
                ksp.solve(v2,v)        # v = (I-0.5*dt*M) \ v2
                t = kt*dt
                if t > nextt:
                    nextt += 0.1
                    PETSc.Sys.Print('time is {0}'.format(t))

        elif cpm == 4:   # SBDF-2: [1-2/3*dt*M]u^{n+1}=4/3*u^n-1/3*u^{n-1}
            dt = 0.1*band.dx
            numtimesteps = int( Tf // dt + 1 )
            v1 = v.copy()
            v1.setArray(uexactfn(dt,band.cp))
            # first do some steps of Forward Euler
#            timesteps_initial = int( 1 // band.dx  )
#            dt_initial = dt / timesteps_initial
#            M1 = M.copy()
#            M1.scale(dt_initial)
#            for k in xrange(timesteps_initial):
#                M1.multAdd(v,v,v1)
#                v1.swap(v)

            # Or do one step of implicit Euler
            # TODO: do we have to sepcify another A and ksp solver?
#            A1 = I.copy()
#            A1.axpy(-dt,M)
#            A1.setUp()
#            ksp1 = PETSc.KSP().create()
#            ksp1.setOperators(A1)
#            ksp1.setType(ksp1.Type.GMRES)
#            pc1 = ksp1.getPC()
#            pc1.setType('none')
#            ksp1.solve(v,v1)
            
            A = I.copy()
            #A.axpy(-2./3*dt,M)
            A = I - 2./3*dt*M
            A.setUp()
            ksp = PETSc.KSP().create()
            ksp.setOperators(A)
            ksp.setType(ksp.Type.GMRES)
            pc = ksp.getPC()
            # slow
            #pc.setType(pc.Type.GAMG)
            # default preconditioner, fastest of all , and seems to be enough accurate...
            pc.setType('none')
            for kt in xrange(2,numtimesteps):
                v2 = 4./3*v1 - v/3
                ksp.solve(v2,v2)        
                v = v1
                v1 = v2
                t = kt*dt
                if t > nextt:
                    nextt += 0.1
                    PETSc.Sys.Print('time is {0}'.format(t))

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
            print('maximal is {0}'.format(ee))
        PETSc.Sys.Print('==================================')   
        del band,v,Eplot,cv 
    if comm.rank == 0:
        import pickle
        ferror = open('error.pickle','w')
        fdx = open('dx.pickle','w')
        pickle.dump(error, ferror)
        pickle.dump(dx,fdx)
    
    
    
