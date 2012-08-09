'''
Created on Aug 2, 2012

@author: nullas
'''
from __future__ import division
import unittest
from numpy import array as a
import numpy.testing as npt
from surfaces import Sphere


class Test(unittest.TestCase):


    def setUp(self):
        
        self.surface = Sphere()


    def tearDown(self):
        pass


    def testSurface(self):
        '''test Surfaces.Sphere'''
        surface = self.surface
        example = ((a([0,0,0]),a([1,0,0]),1),
                   ( a([[-1,0,0],[0,-1,0]]),a([[-1,0,0],[0,-1,0]]),a([0,0]) ),
                   ( a([[0,0,0]]),a([[1,0,0]]),a([1]) ),
                   ( a([1.1,0,0]),a([1,0,0]),0.1))
        for x,y,rd in example:
            rlt,d,_,_ = surface.cp(x)
            npt.assert_array_equal(y,rlt)
            npt.assert_array_equal(d,rd)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()