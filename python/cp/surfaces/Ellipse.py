"""
Closest point representation of an ellipse in 2d

Ellipse(center, a, b)
Default centered at [0.0 0.0] with major axis a = 1.5 and minor axis b = 0.75

Internally, uses cpParamCurveClosed with Newton solves to find cp

"""

import numpy as np

from ParamCurveClosed import ParamCurveClosed


class Ellipse(ParamCurveClosed):
    def __init__(self, center=np.array([0.0, 0.0]), a=1.5, b=0.75):
        self.center = center
        self.a = a
        self.b = b

        # parameterised curve:
        xs = lambda t: self.a * np.cos(t) + self.center[0]
        ys = lambda t: self.b * np.sin(t) + self.center[1]

        # derivative of parametrisation:
        xp = lambda t: -self.a * np.sin(t);
        yp = lambda t: self.b * np.cos(t);

        # second derivative:
        xpp = lambda t: -self.a * np.cos(t);
        ypp = lambda t: -self.b * np.sin(t);

        # start and endpts of parameter variable
        endpt1 = 0.0
        endpt2 = 2*np.pi

        ParamCurveClosed.__init__(self, xs, ys, xp, yp, xpp, ypp, endpt1, endpt2)
