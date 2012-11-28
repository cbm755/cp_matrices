"""
Closest point representation of an ellipsoid in 3d

Ellipsoid(center, a, b)
center = [xc, yc, 0]
Default centered at [0.0 0.0 0.0], rotated around the 'x' axis with major axis a = 1.5 and minor axis b = 0.75

# TODO implement shift center in z-direction (for now set to 0)

"""

import numpy as np

from cp.surfaces.SurfOfRevolution import SurfOfRevolution
from cp.surfaces.Ellipse import Ellipse

class Ellipsoid(SurfOfRevolution):
    def __init__(self, center=np.array([0.0, 0.0, 0.0]), a=1.5, b=0.75, whichaxis='x'):
        self.center = center
        self.a = a
        self.b = b
        # shift to origin, rotate about axis
        SurfOfRevolution.__init__(self, Ellipse, whichaxis, 1, np.array([0.0, 0.0, 0.0]), a, b)
        # shift back
        self.bounding_box = self.bounding_box + self.center
        
    def closest_point(self, xx):
        cpxx, dist, _, _ = SurfOfRevolution.closest_point(self, xx-self.center)
        return cpxx, dist, np.zeros(xx.ndim), {}

    def parametric_grid(self, rez=64):
        xx, yy, zz = SurfOfRevolution.parametric_grid(self,rez)
        xx = xx + self.center[0]
        yy = yy + self.center[1]
        zz = zz + self.center[2]
        return xx, yy, zz
                        
        

