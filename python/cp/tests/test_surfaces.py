import numpy as np
import numpy.testing as npt

from cp.surfaces.Hemisphere import Hemisphere


class TestHemisphere:
    def test_cp_2d(self):
        pass
    def test_cp_3d(self):
        points = np.array([[0, 0, 0],
                           [0, 0, -2],
                           [0, 0, 2],
                           [0.5, 0.5, -1],
                           [2, 2, 0]])
        h = Hemisphere()
        cp, dist, bdy, _ = h.closestPointToCartesianVec(points)
        for cp_i, dist_i, bdy_i, points_i in zip(cp, dist, bdy, points):
            cp_i2, dist_i2, bdy_i2, _ = h.closestPointToCartesian(points_i)
            for x, y in zip((cp_i, dist_i, bdy_i), (cp_i2, dist_i2, bdy_i2)):
                npt.assert_allclose(x, y)
