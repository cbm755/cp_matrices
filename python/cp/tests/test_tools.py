import numpy as np
import numpy.testing as npt

import os

class test_load_ply:
    from cp.tools.io import load_ply
    v, f = load_ply(os.path.join(os.path.split(os.path.abspath(__file__))[0],
                                 "data",
                                 "Armadillo_ascii_converted_using_meshlab.ply"))
    npt.assert_equal(f.max(), len(v) - 1)
    npt.assert_equal(np.setdiff1d(f.ravel(), np.arange(len(v))), np.array([]))

    v2, f2 = load_ply(os.path.join(os.path.split(os.path.abspath(__file__))[0],
                                   "data",
                                   "Armadillo.ply"))
    # big absolute tolerance because meshlab lost many digits when
    # converting to ascii
    npt.assert_allclose(v2, v, atol=1e-4)
    npt.assert_allclose(f2, f)
