"""
Closest point function for a sphere.

Works in multidimensions.
"""
import numpy as np

from Surface import Surface


class Sphere(Surface):
    def __init__(self, center=np.array([0.0, 0.0, 0.0]), radius=1.0):
        self.center = center
        self.radius = radius
        self.dim = len(center)
        self.bounding_box = [center - radius, center + radius]
        if (self.dim == 2) or (self.dim == 3):
            self._hasParam = True

    def closest_point(self, x):
        """TODO: The actual implementation has not been cleaned up.

        It has to be a function that accepts arrays of shape
        (n_points, dim)
        """
        return self.closestPointToCartesian(x)

    def closestPointToCartesian(self, x):
        """
        This implementation works in general dimensions.  But not very
        "clean" and one has to be careful with the details (copy is
        necessary for the computation of distance for example)
        """
        SURF_CEN = self.center
        SURF_SR = self.radius
        
        if x.ndim == 1:
            cpx = x-SURF_CEN
            r = np.linalg.norm(cpx)
            if r == 0:
                cpx[0] += 1
                r = 1
            cpx = SURF_SR/r*cpx + SURF_CEN
            d = np.linalg.norm(cpx-x)
            return cpx,d,0,{}
            
        cpx = x-SURF_CEN
        r = np.sqrt(np.sum(pow(cpx,2),axis=1))
        if (r.any()==0):
            ind = np.where(r == 0)
            cpx[ind,0] = SURF_SR
            r[ind] = SURF_SR


        r = SURF_SR / r
        for dim in range(x.shape[1]):
            cpx[:,dim] *= r
        cpx += SURF_CEN
        dist = np.sqrt(np.sum(pow(cpx - x, 2),axis=1))
        return cpx, dist, np.zeros(x.ndim), {}

    def parametric_grid(self, rez=20):
        """ Return a mesh, for example for plotting with mlab
        """
        # parametric variables
        #u=np.r_[0:2*pi:10j]
        #v=np.r_[0:0.5*pi:10j]
        th = np.linspace(0, 2*np.pi, 2*rez, True)
        phi = np.linspace(0, 1*np.pi, rez, True)
        x = self.radius*np.outer(np.cos(th), np.sin(phi)) + self.center[0]
        y = self.radius*np.outer(np.sin(th), np.sin(phi)) + self.center[1]
        z = self.radius*np.outer(np.ones(th.size), np.cos(phi)) + self.center[2]
        # TODO: return a list?  In case we have multiple components
        return x, y, z
