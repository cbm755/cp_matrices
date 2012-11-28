"""
Testing ParamCurveClosed
Build a parameterised circle, and check that it gives the same results as the Circle code, for grid, distance, cp
"""

import numpy as np
from cp.surfaces import Circle
from cp.surfaces.ParamCurveClosed import ParamCurveClosed

import matplotlib.pyplot as plt

# Using circle function

s1 = Circle()

p = 3
diff_stencil_arm = 1
dim = 2

cp1, distance1, grid1, dx1 = s1.grid(num_blocks_per_dim=41,
                                   levels=1,
                                   p=p,
                                   diff_stencil_arm=diff_stencil_arm)


# Parameterised Circle - using ParamCurveClosed


xs = lambda t: np.cos(t)
ys = lambda t: np.sin(t)
xp = lambda t: -np.sin(t)
yp = lambda t: np.cos(t)
xpp = lambda t: -np.cos(t)
ypp = lambda t: -np.sin(t)

endpt1 = 0.0
endpt2 = 2*np.pi

s2 = ParamCurveClosed(xs,ys,xp,yp,xpp,ypp,endpt1,endpt2)

xp, yp = s2.parametric_grid(65)
plt.plot(xp,yp,'x')
plt.axis('equal')
#plt.show()


cp2, distance2, grid2, dx2 = s2.grid(num_blocks_per_dim=41,
                                   levels=1,
                                   p=p,

                                   diff_stencil_arm=diff_stencil_arm)



# Check that the virtual grids are the same

# Circle virtual grid
ll1 = np.array(dim * [grid1.min()]) - 3 * dx1
ur1 = np.array(dim * [grid1.max()]) + 3 * dx1
virtual_grid_shape1 = np.abs(ur1-ll1) / dx1 + 1
print "Grids: lower left, upper right, shape"
print "Circle:", ll1, ur1, virtual_grid_shape1

# ParamCurve virtual grid
ll2 = np.array(dim * [grid2.min()]) - 3 * dx2
ur2 = np.array(dim * [grid2.max()]) + 3 * dx2
virtual_grid_shape2 = np.abs(ur2-ll2) / dx2 + 1
print "ParamCurve:", ll2, ur2, virtual_grid_shape2

diffgrid = abs(grid2 - grid1)
print "grid difference", diffgrid.max()

# Check that the closest points are the same 

diffcp = abs(cp2 - cp1)
print "cp difference", diffcp.max()

diffdist = abs(distance2 - distance1)
print "distance difference", diffdist.max()
