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
def cpz(cp):
    return cp[:,2]

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
    
def monitor(ts, i, t, x):
    PETSc.Sys.Print(
        'At {0} max-norm is {1}'.format(t, x.norm(PETSc.NormType.INFINITY)))

if __name__ == '__main__':
    MBlocklist = [20,40,80,160]
    error = []
    dx = []
    ft = []
    l = sp.linspace(-sp.pi, sp.pi, 100)
    x,y = sp.mgrid[-sp.pi:sp.pi:100j,-sp.pi/2:sp.pi/2:50j]
    x = x.reshape((-1,))
    y = y.reshape((-1,))
    points = sp.column_stack((sp.cos(x)*sp.cos(y),sp.sin(x)*sp.cos(y),sp.sin(y)))
    exactu = sp.sin(y)
    
    comm = MPI.COMM_WORLD
    
    
    vsize = points.shape[0]
    vAssigned = vsize // comm.size + int(comm.rank < (vsize % comm.size))


    vstart = comm.exscan(vAssigned)
    if comm.rank == 0:
        vstart = 0
    for MBlock in MBlocklist:
        opt = {'M':MBlock,'m':4,'d':3}
        surface = Sphere()
        
        band = Band(surface,comm,opt)
        la,lv,gv,wv = band.createGLVectors()
        v = band.getCoordinates() 
        vv = sp.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
        weights = sp.array([1,1,1,1,1,1])*(1/band.dx**2)
        L = band.createAnyMat(vv, weights, (7,3))
        PETSc.Sys.Print('Laplacian')

        M = band.createExtensionMat()
        PETSc.Sys.Print('ExtensionMat built')

        ts = PETSc.TS().create(comm=comm)
        ts.setProblemType(ts.ProblemType.LINEAR)
        ts.setType(ts.Type.BEULER)
        


        PETSc.Sys.Print('ML creating')

        LM = M.matMult(L)

        I = PETSc.Mat().create(comm=comm)
        I.setSizes((gv.sizes,gv.sizes))
        I.setFromOptions()
        I.setPreallocationNNZ(1)
        wv.set(1)
        I.setDiagonal(wv)
        I.assemble()

        
        LM.axpy(-6 / band.dx ** 2, I)
        LM.mult(wv,gv)
        
        PETSc.Sys.Print('ML created')

        ts.setTime(0.0)
        ts.setTimeStep(0.01 * band.dx)
        ts.setDuration(1)
        gv.set(1)
        
        LM.mult(gv,wv)

        
        PETSc.Sys.Print('wv inf-norm is {0}'.format(wv.norm(PETSc.NormType.INFINITY)))
        
        
        band.initialu(cpz)
        ts.setSolution(gv)
        ts.setFromOptions()
        ts.setRHSFunction(PETSc.TS.computeRHSFunctionLinear,wv)
        ts.setRHSJacobian(PETSc.TS.computeRHSJacobianConstant,LM,LM)
        ts.setMonitor(monitor)
#        ts.step()
        PETSc.Sys.Print('Begin to solve')
        t = ts.solve(gv)
             

        mv = band.createExtensionMatForLoop(cp=points[vstart:vstart+vAssigned])
        cv = mv.getVecLeft()
        mv.mult(gv,cv)
        cv = band.toZeroStatic(cv)
        if comm.rank == 0:
            cv = cv.getArray()
            exu = sp.exp(-2*t)*exactu
            ee = norm(cv-exu, sp.inf)
            error.append(ee)
            dx.append( band.dx )
            
        M.mult(gv,wv)
        
    
        
        PETSc.Sys.Print('==================================')  
        PETSc.Sys.Print('wv inf-norm is {0}'.format(wv.norm(PETSc.NormType.INFINITY))) 
        if comm.rank == 0: 
            print('maximal is {0}'.format(ee))
            
        PETSc.Sys.Print('==================================')   
        del band,v,mv,cv,M,L,LM,ts 
    if comm.rank == 0:
        import pickle
        ferror = open('error.pickle','w')
        fdx = open('dx.pickle','w')
        pickle.dump(error, ferror)
        pickle.dump(dx,fdx)
    
    
    