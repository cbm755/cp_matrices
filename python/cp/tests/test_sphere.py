'''
Created on Aug 2, 2012

@author: nullas
'''
from __future__ import division
import unittest
from numpy import array as a
import numpy.testing as npt
import scipy as sp
from cp.surfaces import Sphere


class Test(unittest.TestCase):


    def setUp(self):
        
        self.surface = Sphere()


    def tearDown(self):
        pass


    def testSurface(self):
        '''test Surfaces.Sphere'''
        surface = self.surface
        example = (( a([0,0,0]),a([1,0,0]),1 ),
                   ( a([2,2,2]),a([1,1,1])/sp.sqrt(3),2*sp.sqrt(3)-1),
                   ( a([[-1,0,0],[0,-1,0]]),a([[-1,0,0],[0,-1,0]]),a([0,0]) ),
                   ( a([[0,0,0]]),a([[1,0,0]]),a([1]) ),
                   ( a([1.1,0,0]),a([1,0,0]),0.1))
        for x,y,rd in example:
            rlt,d,_,_ = surface.cp(x)
            npt.assert_allclose(y,rlt)
            npt.assert_allclose(d,rd)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()