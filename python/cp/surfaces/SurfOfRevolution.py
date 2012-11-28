"""
Surface of Revolution, a curve rotated to make a surface
Implemented using a given curve (Circle, Ellipse, ParamCurveClosed etc), rotated around the x or y axis

Default
Ellipse, rotated around the x axis

Input
- Surface class
- whichaxis, axis rotated about (x or y)
- symm, True if object is symmetric about the axis (eg ellipse), False if defined only in the upper half plane
- cpfdata, optional parameters for the surface class

Disclaimers / warnings apply, restrictions exist on the curve
eg. curve must either be symmetric about the axis or defined in the upper half plane only

TODO: boundaries

"""

import numpy as np
from Surface import Surface
from cp.surfaces.Ellipse import Ellipse
from cp.surfaces.coordinate_transform import cart2pol, pol2cart


class SurfOfRevolution(Surface):
    def __init__(self, cpf=Ellipse, whichaxis='x', symm=True, *cpfdata):
        self.dim = 3
        self.axis = whichaxis
        self.symm = symm
        self.curve = cpf(*cpfdata)
        self.bounding_box = self.find_bounding_box()
        self._hasParam = True



    def find_bounding_box(self):
        ll, ur = self.curve.bounding_box;
        if self.axis == 'x':
            # rotated about x axis
            xmin = ll[0]
            xmax = ur[0]
            ymin = min(-ur[1], ll[1])
            ymax = max(ur[1], -ll[1])
            zmin = ymin
            zmax = ymax
        elif self.axis == 'y':
            # rotated about y axis
            xmin = min(-ur[0], ll[0])
            xmax = max(ur[0], -ll[0])
            ymin = ll[1]
            ymax = ur[1]
            zmin = xmin
            zmax = xmax
        else:
            print "Error: axis of revolution should be x or y"

        return [np.array([xmin, ymin, zmin]), np.array([xmax, ymax, zmax])]


    def closest_point(self, xx):
        x,y,z = xx.T
        # no boundaries for now
        if self.axis == 'x':
            th, r, zz = cart2pol(y,z,x)
            cpzzr, dist, _, _ = self.curve.closest_point(np.column_stack((zz, r)))
            cpzz, cpr = cpzzr.T
            cpy, cpz, cpx = pol2cart(th, cpr, cpzz)
        elif self.axis == 'y':
            th, r, zz = cart2pol(z,x,y)
            cprzz, dist, _, _ = self.curve.closest_point(np.column_stack((r, zz)))
            cpr, cpzz = cprzz.T
            cpz, cpx, cpy = pol2cart(th, cpr, cpzz)
        
        cpxx = np.column_stack((cpx, cpy, cpz))
        return cpxx, dist, np.zeros(xx.ndim), {}                                       
        


    def parametric_grid(self, rez=64):
        """ Parameterized form (for plotting with surf) """
        # TODO: this depends on whether the original curve was defined in the whole or half plane
        # Should probably be overwritten by the particular implementation of the SurfOfRevolution for best results 
        X, Y = self.curve.parametric_grid(rez)
        if self.symm:
            # defined in whole plane
            uendpt = np.pi
        else:
            # defined in half plane
            uendpt = 2*np.pi
        u = np.linspace(0.0, uendpt, rez)

        if self.axis == 'x':
            xx = np.tile(X,(rez,1)).T
            yy = np.outer(Y, np.cos(u))
            zz = np.outer(Y, np.sin(u))
        elif self.axis == 'y':
            xx = np.outer(X, np.cos(u))
            yy = np.tile(Y,(rez,1)).T
            zz = np.outer(X, np.sin(u))
        else:
            print "Error: axis of revolution should be x or y"

        return xx, yy, zz
 
