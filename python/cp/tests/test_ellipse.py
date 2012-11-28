"""
Testing Ellipse
"""

import numpy as np
from cp.surfaces.Ellipse import Ellipse

import matplotlib.pyplot as plt


s = Ellipse()

p = 3
diff_stencil_arm = 1
dim = 2

cp, distance, grid, dx = s.grid(num_blocks_per_dim=41,
                                   levels=1,
                                   p=p,
                                   diff_stencil_arm=diff_stencil_arm)


xp, yp = s.parametric_grid(65)
plt.plot(xp,yp,'b')
plt.plot(cp[:,0],cp[:,1],'gx')
plt.axis('equal')
plt.show()

# virtual grid
ll = np.array(dim * [grid.min()]) - 3 * dx
ur = np.array(dim * [grid.max()]) + 3 * dx
virtual_grid_shape = np.abs(ur-ll) / dx + 1
