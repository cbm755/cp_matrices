import numpy as np
import numpy.testing as npt

import os

def test_load_ply():
    from cp.tools.io import load_ply
    v, f = load_ply(os.path.join(os.path.split(os.path.abspath(__file__))[0],
                                 "data",
                                 "Armadillo_ascii_converted_using_meshlab.ply"))
    npt.assert_equal(f.max(), len(v) - 1)
    npt.assert_equal(np.setdiff1d(f.ravel(), np.arange(len(v))), np.array([]))
