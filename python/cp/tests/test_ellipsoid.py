"""
Testing Ellipse
"""

import numpy as np
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab

from cp.surfaces.Ellipsoid import Ellipsoid




s = Ellipsoid(np.array([3.0, 4.0, 0.0]), 2., 3.)

p = 3
diff_stencil_arm = 1
dim = 2

cp, distance, grid, dx = s.grid(num_blocks_per_dim=41,
                                   levels=1,
                                   p=p,
                                   diff_stencil_arm=diff_stencil_arm)

xx, yy, zz = s.parametric_grid(65)

mlab.mesh(xx,yy,zz)
mlab.show()


