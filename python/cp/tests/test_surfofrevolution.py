"""
Testing Surface of Revolution
Build a Circle of radius 2 rotated around the x axis, and compare to Sphere
"""

import numpy as np
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab


from cp.surfaces import Sphere
from cp.surfaces import Circle
from cp.surfaces.SurfOfRevolution import SurfOfRevolution


# NB must import the surface you want to use

p = 3
diff_stencil_arm = 1
dim = 3

# parameterised sphere
s1 = SurfOfRevolution(Circle, 'x', 1, np.array([0.0, 0.0]), 2)


cp1, distance1, grid1, dx1 = s1.grid(num_blocks_per_dim=41,
                                   levels=1,
                                   p=p,
                                   diff_stencil_arm=diff_stencil_arm)
# Sphere
s2 = Sphere(np.array([0.0, 0.0, 0.0]), 2)


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



# Plot surface
xx, yy, zz = s1.parametric_grid()

mlab.mesh(xx,yy,zz)
mlab.show()


# plot a torus
#s3 = SurfOfRevolution(Circle, 'x', 0, np.array([3, 3]))

#xx, yy, zz = s3.parametric_grid()

#mlab.mesh(xx,yy,zz)
#mlab.show()
