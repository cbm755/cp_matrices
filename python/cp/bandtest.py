'''
Created on Jul 30, 2012

@author: nullas
'''
from __future__ import division

import unittest
import exceptions

from surfaces import Sphere as Surface
from petsc.band import Band
from mpi4py import MPI


try:
    import petsc4py
    import sys
    from petsc4py import PETSc
    petsc4py.init(sys.argv)
    from numpy import array as a
    import numpy.testing as npt
except Exception as exp:
    print exp

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
        example = ( ( 0, a([-1.9,-1.9,-1.9])), 
                    ( 1, a([-1.7,-1.9,-1.9])) )
        for ind,right in example:
            sub = self.bnd.BlockInd2SubWithoutBand(ind)
            rslt = self.bnd.BlockSub2CenterCarWithoutBand(sub)
            npt.assert_array_equal(rslt, right)

            
    def testVisually(self):
        '''visually blocks selected.'''
#        if self.comm.rank == 0:
        g2z,zvec = PETSc.Scatter().toZero(self.bnd.gsubBlockWBand)
        g2z.scatter(self.bnd.gsubBlockWBand,zvec, PETSc.InsertMode.INSERT)
        x = self.bnd.BlockSub2CenterCarWithoutBand(zvec.getArray())
        lx = self.bnd.BlockSub2CenterCarWithoutBand(self.bnd.gsubBlockWBand.getArray())
        print x
        import pylab as pl
        from mpl_toolkits.mplot3d import Axes3D #@UnusedImport
        fig = pl.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter3D(x[::self.bnd.Dim],x[1::self.bnd.Dim],x[2::self.bnd.Dim],c='blue',marker='o')
        ax.scatter3D(lx[::self.bnd.Dim],lx[1::self.bnd.Dim],lx[2::self.bnd.Dim],c='red',marker='D')
        pl.savefig('testVis.png')
        
            


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()