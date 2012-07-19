import numpy as np
import numpy.testing

def test_LagrangeWeights1D():
    from cp import cpOps
    numpy.testing.assert_equal(cpOps.LagrangeWeights1D(-0.25, np.array([0.]), 0.25, 4), np.array([[0, 1, 0, 0]]))
