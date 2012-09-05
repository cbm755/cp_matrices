"""
Closest point representation of a circle.

A special case of Sphere but implemented differently because this
version is clearer to understand.
"""
import numpy as np

from Surface import Surface
from coordinate_transform import cart2pol, pol2cart


class Circle(Surface):
    def __init__(self, center=np.array([0.0, 0.0]), radius=1.0):
        # TODO: could make it subclass of sphere and just call:
        # self.super(self, center=center, radius=radius)
        self.center = center
        self.radius = radius
        self.dim = 2
        self.bounding_box = [center - radius, center + radius]
        self._hasParam = True

    def closest_point(self, x):
        th, r = cart2pol(*(x-self.center).T)
        cpx, cpy = pol2cart(th, self.radius)
        dist = np.abs(r - self.radius)
        cp = np.column_stack((cpx, cpy)) + self.center
        # Maybe we should return np.zeros(dist.shape) instead of just one 0?
        return cp, dist, 0, {}
    
    def parametric_grid(self, rez=256):
        """
        Parameterized form (for plotting)
        """
        th = np.linspace(0, 2*np.pi, num=rez)
        circ = self.radius * np.exp(1j*th)
        X = np.real(circ) + self.center[0]
        Y = np.imag(circ) + self.center[1]
        return X, Y
