'''
Created on Aug 10, 2012

@author: nullas
'''
from __future__ import division
import scipy as sp
from cp.tools.io import load_ply
from cp.mesh_cp import MeshCP

class MeshWrapper(object):
    '''
    Wrapper for handling mesh.
    '''


    def __init__(self,ff = None):
        '''
        Constructor
        '''
        if ff is None:
            ff = open('Armadillo_ascii_converted_using_meshlab.ply')
        v,f = load_ply(ff)
        ma = v.max()
        mi = v.min()
        med = (ma+mi)/2
        v -= med
        v /= (ma-mi)/2
        
        
        self.v = v
        self.f = f
        self.meshcp = MeshCP(v,f)
        
    def cp(self,points,bounds = None):
        if bounds is None:
            return self.meshcp.query(points),0,0,0
        else:
            return self.meshcp.query(points, bounds),0,0,0
        
        
        