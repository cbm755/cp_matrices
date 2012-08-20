import numpy as np
import numpy.testing as npt

from cp.surfaces.Hemisphere import Hemisphere


class TestHemisphere:
    def test_cp_2d(self):
        h = Hemisphere(np.zeros(2))
        points = np.array([[0, 0],
                           [-1, -1],
                           [-1, 1],
                           [2, 2]])
        cp, dist, bdy, _ = h.closestPointToCartesianVec(points)
        for cp_i, dist_i, bdy_i, points_i in zip(cp, dist, bdy, points):
            cp_i2, dist_i2, bdy_i2, _ = h.closestPointToCartesian(points_i)
            self.check_cp_results((cp_i, dist_i, bdy_i),
                                  (cp_i2, dist_i2, bdy_i2))

    def check_cp_results(self, res_a, res_b):
        for x, y in zip(res_a, res_b):
            npt.assert_allclose(x, y)

    def test_cp_3d(self):
        points = np.array([[0, 0, 0],
                           [0, 0, -2],
                           [0, 0, 2],
                           [0.5, 0.5, -1],
                           [2, 2, 0]])
        h = Hemisphere()
        cp, dist, bdy, _ = h.closestPointToCartesianVec(points)
        for cp_i, dist_i, bdy_i, points_i in zip(cp, dist, bdy, points):
            # If only we were using py3... "*res_a, _ = "
            cp_i2, dist_i2, bdy_i2, _ = h.closestPointToCartesian(points_i)
            self.check_cp_results((cp_i, dist_i, bdy_i),
                                  (cp_i2, dist_i2, bdy_i2))
