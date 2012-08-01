import numpy as np
import numpy.testing as npt

from cp.cpOps import LagrangeWeights1D


class test_LagrangeWeights1D():
    def test_nodes(self):
        n = 5
        dx = 0.25
        for i in xrange(n):
            expected = np.zeros(n)
            expected[i] = 1
            result = LagrangeWeights1D(0, i * dx, dx, n)
            yield self.check_result, expected, result
            
    def check_result(self, expected, result):
        npt.assert_allclose(expected, result)

    def test_many_points(self):
        npt.assert_allclose(LagrangeWeights1D(0, 0.1 * np.arange(4), 0.1, 5),
                            np.eye(4,5))

    def test_many_points_and_basepoints(self):
        result = np.zeros((4,5))
        result[:, 0] = 1
        npt.assert_allclose(LagrangeWeights1D(0.1 * np.arange(4),
                                          0.1 * np.arange(4), 0.1, 5),
                            result)

    def test_many_points_and_basepoints_and_dx(self):
        x = np.array([[0, 1, 2],
                      [0, 0, 0]])
        xg = np.array([[0.1,  1.2,  2.3],
                      [0.2,  1. ,  0.6]])
        dx = np.array([0.1, 0.2, 0.3])
        npt.assert_allclose(LagrangeWeights1D(x, xg, dx, 6),
                            np.array([[[-0.,  1.,  0., -0.,  0., -0.],
                                       [-0.,  1.,  0., -0.,  0., -0.],
                                       [-0.,  1.,  0., -0.,  0., -0.]],
                                       
                                      [[ 0., -0.,  1.,  0., -0.,  0.],
                                       [-0.,  0., -0.,  0., -0.,  1.],
                                       [ 0., -0.,  1.,  0., -0.,  0.]]]))
