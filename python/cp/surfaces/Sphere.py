"""
Closest point function for a sphere.

Works in multidimensions.
"""
import numpy as np

from Surface import Surface
from coordinate_transform import sph2cart

class Sphere(Surface):
    def __init__(self, center=np.array([0.0, 0.0, 0.0]), radius=1.0):
        self._center = center
        self._radius = radius
        self._dim = len(center)
        self._bb = [center - radius, center + radius]
        if self._dim == 2 or self._dim == 3:
            self._hasParam = True
        # TODO: superclass knows about dimension?
        #super(Hemisphere, self).__init__()

    def closestPointToCartesian(self, x):
        """
        This implementation works in general dimensions.  But not very
        "clean" and one has to be careful with the details (copy is
        necessary for the computation of distance for example)
        """
        SURF_CEN = self._center
        SURF_SR = self._radius

        r = np.linalg.norm(x - SURF_CEN, 2)
        if r==0:
            tx = SURF_CEN[0] + 1.0
            r = 1.0
        else:
            tx = x[0]
        xm = x.copy()
        xm[0] = tx

        c = SURF_SR / r
        cpx = c*(xm-SURF_CEN) + SURF_CEN
        dist = np.linalg.norm(cpx - x, 2)
        return cpx, dist, 0, {}

    cp = closestPointToCartesian

    def ParamGrid(self, rez=20):
        """ Return a mesh, for example for plotting with mlab
        """
        rad = self._radius
        cen = self._center
        # parametric variables
        #u=np.r_[0:2*pi:10j]
        #v=np.r_[0:pi:10j]
        th = np.linspace(0, 2*np.pi, 2*rez)
        # Since sph2cart uses the Matlab convention, elevation angle
        # goes from -pi/2 to pi/2
        phi = np.linspace(-np.pi/2, np.pi/2, rez)
        TH, PHI = np.meshgrid(th, phi)
        x, y, z = sph2cart(TH, PHI, rad)
        # TODO: return a list?  In case we have multiple components
        return x + cen[0],y + cen[1], z + cen[2]
