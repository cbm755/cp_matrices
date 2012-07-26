import numpy as np
import numpy.testing as npt

def test_LagrangeWeights1D():
    from cp.cpOps import LagrangeWeights1D
    npt.assert_allclose(LagrangeWeights1D(-0.25, 0., 0.25, 4),
                        np.array([0, 1, 0, 0]))
    npt.assert_allclose(LagrangeWeights1D(0, 0.2, 0.1, 4),
                        np.array([ 0., -0.,  1.,  0.]))
    npt.assert_allclose(LagrangeWeights1D(0, 0.1 * np.arange(4), 0.1, 4),
                        np.array([[ 1.,  0., -0.,  0.],
                                  [-0.,  1.,  0., -0.],
                                  [ 0., -0.,  1.,  0.],
                                  [-0.,  0., -0.,  1.]]))
    npt.assert_allclose(LagrangeWeights1D(0.1 * np.arange(4),
                                          0.1 * np.arange(4), 0.1, 4),
                        np.array([[ 1.,  0., -0.,  0.],
                                  [ 1.,  0., -0.,  0.],
                                  [ 1.,  0., -0.,  0.],
                                  [ 1.,  0., -0.,  0.]]))
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
