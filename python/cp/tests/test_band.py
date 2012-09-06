'''
Created on Jul 30, 2012

@author: nullas
'''
from __future__ import division

import unittest
import exceptions
from nose.plugins.attrib import attr

from cp.surfaces import Sphere as Surface
from cp.petsc.band import Band
from mpi4py import MPI
import petsc4py
import sys
from petsc4py import PETSc
petsc4py.init(sys.argv)
from numpy import array as a
import numpy.testing as npt


class TestBand(unittest.TestCase):


    def setUp(self):
        self.comm = MPI.COMM_WORLD
        surface = Surface()
        self.bnd = Band(surface,self.comm)
        self.bnd.createGLVectors()

    def tearDown(self):
        pass


    def testBlockInd2SubWithoutBand(self):
        '''test BlockInd2SubWithoutBand'''
        if 1 == 1:return 
        example = ( (5, a([5,0,0]) ),
                    (21,a([1,1,0]) ),
                    (a([401,402]),a([[1,0,1],[2,0,1]])) )
        for ind,sub in example:
            result = self.bnd.BlockInd2SubWithoutBand(ind)
            npt.assert_array_equal(result,sub,'{0} is wrong '.format(ind))
        ExceptionExample = (10000,8000)
        for i in ExceptionExample:
            self.assertRaises(exceptions.IndexError, self.bnd.BlockInd2SubWithoutBand,i)
            
    def testBlockSub2CenterCarWithoutBand(self):
        '''test BlockSub2CenterCarWithoutBand'''
        if 1 == 1: return
        example = ( ( 0, a([-1.9,-1.9,-1.9])), 
                    ( 1, a([-1.7,-1.9,-1.9])) )
        for ind,right in example:
            sub = self.bnd.BlockInd2SubWithoutBand(ind)
            rslt = self.bnd.BlockSub2CenterCarWithoutBand(sub)
            npt.assert_array_equal(rslt, right)

   
    @attr(visual=1)
    def testVisually(self):
        '''blocks selected visually.'''
#        if self.comm.rank == 0:
        g2z,zvec = PETSc.Scatter().toZero(self.bnd.gindBlockWBand)
        g2z.scatter(self.bnd.gindBlockWBand,zvec, PETSc.InsertMode.INSERT)
        x = self.bnd.BlockSub2CenterCarWithoutBand(\
                                                   self.bnd.BlockInd2SubWithoutBand(zvec.getArray()) )
        lx = self.bnd.BlockSub2CenterCarWithoutBand(\
                                                    self.bnd.BlockInd2SubWithoutBand(self.bnd.gindBlockWBand.getArray()))
        try:
            try:
                from mayavi import mlab
            except ImportError:
                from enthought.mayavi import mlab

            if self.comm.rank == 0:
                mlab.figure()
                mlab.points3d(x[:,0],x[:,1],x[:,2])
            mlab.figure()
            mlab.points3d(lx[:,0],lx[:,1],lx[:,2])
            mlab.show()
            #fig.add(pts1)
            #fig.add(pts2)
        except ImportError:
            import pylab as pl
            from mpl_toolkits.mplot3d import Axes3D #@UnusedImport
            fig = pl.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter3D(x[:,0],x[:,1],x[:,2],c='blue',marker='o')
            ax.scatter3D(lx[:,0],lx[:,1],lx[:,2],c='red',marker='D')
            pl.savefig('testVis{0}.png'.format(self.comm.rank))
            pl.show()
        
    def testNorm1(self):
        '''test Norm1.'''
        example = ( ( a( [1,1,1] ), 1 ) ,
                    ( a( [[1,3,4],[2,-2,1]] ), a( [4,2] )) 
                   )
        for x,y in example:
            rslt = Band.norm1(x)
            npt.assert_allclose(rslt, y)
            
            
    def testgetCoordinatesWithGhost(self):
        '''What the return values of this function look like.'''
        #This test takes time.
        if 1 == 1 :return
        x = self.bnd.getCoordinatesWithGhost()
        try:
            from mayavi import mlab
            mlab.figure()
            mlab.points3d(x[:,0],x[:,1],x[:,2])
            #fig.add(pts)
            mlab.show()
            
        except ImportError:
            import pylab as pl
            from mpl_toolkits.mplot3d import Axes3D #@UnusedImport
            fig = pl.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter3D(x[:,0],x[:,1],x[:,2])
            pl.savefig('x{0}.png'.format(self.comm.rank))
            pl.show()
            
    def testgetCoordinates(self):
        '''What the return values of this function look like.'''
        #This test takes time.
        if 1 == 1 :return
        x = self.bnd.getCoordinates()
        try:
            from mayavi import mlab
            mlab.figure()
            mlab.points3d(x[:,0],x[:,1],x[:,2])
            #fig.add(pts)
            mlab.show()
            
        except ImportError:
            import pylab as pl
            from mpl_toolkits.mplot3d import Axes3D #@UnusedImport
            fig = pl.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter3D(x[:,0],x[:,1],x[:,2])
            pl.savefig('x{0}.png'.format(self.comm.rank))
            pl.show()
        
    def testBuildExtMat(self):
        self.bnd.createExtensionMat()
        
    def testCreateAnyMat(self):
        v = a([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
        weights = a([-6,1,1,1,1,1,1])
        self.bnd.createAnyMat(v, weights, (7,3))
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
