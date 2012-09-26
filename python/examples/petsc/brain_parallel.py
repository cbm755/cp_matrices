''' Solve PDEs on a triangulated brain, using Xin's parallel code'''

import numpy as np
import timeit
from mpi4py import MPI
from cp.surfaces.MeshWrapper import MeshWrapper
from cp.surfaces.Sphere import Sphere
from cp.petsc.band_with_doc import Band
import sys
import petsc4py
from petsc4py import PETSc
from cp.surfaces.coordinate_transform import cart2sph

PLOT = False
#PLOT = True

if PLOT:
    try:
        from mayavi import mlab
    except ImportError:
        from enthought.mayavi import mlab

petsc4py.init(sys.argv)

def uinitialfn(cp):
    #th,phi,r = cart2sph(cp[:,0],cp[:,1],cp[:,2])
    #return np.cos(phi + np.pi / 2)
    return cp[:,0]**3

if __name__ == '__main__':
    opt = {'M':40,'m':2,'d':3,'xmin':-2,'xmax':2}
    #opt = {'M':10,'m':4,'d':3}
    comm = MPI.COMM_WORLD

    PETSc.Sys.Print('Starting to load ply file and initialize the mesh wrapper...')
    m = MeshWrapper(ff='brain-lh_scale_1.ply')
    PETSc.Sys.Print('ply file loaded') 

    PETSc.Sys.Print('Starting to select coarse blocks...')
    band = Band(m,comm,opt)
    PETSc.Sys.Print('Coarse blocks have been selected')

    v,v2 = band.createGlobalVectors()
    grid = band.getCoordinates() 

    PETSc.Sys.Print('Starting to compute closest points...')
    band.computeCP()
    PETSc.Sys.Print('Closest Points have been computed')

    PETSc.Sys.Print('Starting to build Laplacian Matrix ...')
    L = band.createLaplacianMat()
    PETSc.Sys.Print('Laplacian built')

    PETSc.Sys.Print('Starting to build Extension Matrix ...')
    E = band.createExtensionMat()
    E1 = band.createExtensionMat(p=1)
    PETSc.Sys.Print('ExtensionMat built')
   
    #v.setArray(uinitialfn(band.cp)) 
    v.setArray(uinitialfn(grid)) 
    PETSc.Sys.Print('Initial value has been set')

    if PLOT:
        PETSc.Sys.Print('Starting to build Eplot ...')
        Eplot = band.createExtensionMat(cp=m.vertices)
        PETSc.Sys.Print('Eplot built')
        vplot = Eplot*v

        # Plotting code. Build a pipeline to be able to change the data later.
        #src = mlab.pipeline.grid_source(xp, yp, zp,
        #                                scalars=(Eplot * u).reshape(xp.shape))
        if comm.rank == 0:
            src = mlab.pipeline.triangular_mesh_source(m.vertices[:,0],
                                                       m.vertices[:,1],
                                                       m.vertices[:,2], m.faces,
                                                       scalars=(vplot.getArray()))
            normals = mlab.pipeline.poly_data_normals(src)
            surf = mlab.pipeline.surface(normals)
            mlab.colorbar()
  
    # choose a closest point method
    cpm = 2
    
    # choose timestep
    if cpm == 0:
        dt = 0.2*band.dx**2
    elif cpm == 1:
        dt = 0.2*band.dx**2
    elif cpm == 2:
        dt = 0.5*band.dx
        
    # various CPM algorithms, if cpm>=1: use various MOLs, vGMM
    if cpm >= 1:  
        lamb = 2*band.Dim / band.dx**2
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
        M.axpy(-lamb,IminusE)

    # Some PETSc set up
    if cpm == 0:
        L *= dt
    if cpm == 1:
        M *= dt
    if cpm == 2:
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

    t = 0
    #Tf = 100*dt
    Tf = 1
    numtimesteps = int( Tf // dt + 1 )

    PETSc.Sys.Print('Begin to solve')
  
    start_time = timeit.default_timer()
    for kt in xrange(numtimesteps):
        if cpm == 0: # Ruuth-Merriman
            L.multAdd(v,v,v2)      # v2 = v + (dt*L)*v 
            E.mult(v2,v)           # v = E*v2
        elif cpm == 1:  # explicit Euler, von Glehn--Maerz--Macdonald
            M.multAdd(v, v, v2)    # v2 = v + (dt*M)*v
            # todo: how to just assign data in v2 to v?
            v2.swap(v)
        elif cpm == 2:   # implicit Euler, vGMM
            ksp.solve(v,v2)
            v = v2
            t = kt*dt
        
        t = kt*dt
        if not kt%100 or kt == (numtimesteps-1):
            PETSc.Sys.Print("time: {0:2f}, {1:2f} %".format(t, 100 * float(kt) / numtimesteps))
            if PLOT:
                vplot = Eplot * v
                if comm.rank == 0:
                    src.data.point_data.scalars = vplot.getArray()
                    src.data.point_data.scalars.name = 'scalars'
                    src.data.modified()

    print "Times, rank=", comm.rank, "time=", timeit.default_timer() - start_time

    PETSc.Sys.Print('maximum of the absolute value of the solution is {0}'.format(v.max()[1]))    

del band,grid   
if PLOT:
    del Eplot

