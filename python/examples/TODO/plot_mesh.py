'''
Created on Aug 15, 2012

@author: nullas
'''

from cp.surfaces.MeshWrapper import MeshWrapper
from mayavi import mlab as pl
import scipy as sp
import petsc4py
petsc4py.init()
from petsc4py import PETSc

def initialu(cp):
    rlt = sp.ones(cp.shape[0])
    (ind,) = sp.where(2.3*cp[:,0]+cp[:,2] > 0)
    rlt[ind] = 0
    return rlt

if __name__ == '__main__':
    surface = MeshWrapper('eight.ply')
    pl.figure(bgcolor=(1,1,1),fgcolor=(0.5,0.5,0.5))
    pl.triangular_mesh(surface.v[:,0],surface.v[:,1],surface.v[:,2],surface.f,scalars=initialu(surface.v))
#    pl.points3d(surface.v[:,0],surface.v[:,1],surface.v[:,2],initialu(surface.v),mode='point')
    pl.colorbar()
    vv=PETSc.Viewer().createBinary('gv.dat','r')
    cv = PETSc.Vec().load(vv)
    cv = cv.getArray()
    pl.figure(bgcolor=(1,1,1),fgcolor=(0.5,0.5,0.5))
    pl.triangular_mesh(surface.v[:,0],surface.v[:,1],surface.v[:,2],surface.f,scalars=cv)
    pl.colorbar()
    pl.show()