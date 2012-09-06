"""Solves the heat equation on a true sphere."""
import numpy as np
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab

from cp.surfaces import Sphere
from cp.build_matrices import build_interp_matrix, build_diff_matrix
from cp.surfaces.coordinate_transform import cart2sph


PLOT = False

s = Sphere()

p = 3
diff_stencil_arm = 1
dim = 3

# As a byproduct of finding the banded grid, we already have its
# closest points, so we don't really have to call s.closest_point()
cp, distance, grid, dx = s.grid(num_blocks_per_dim=41,
                                   levels=1,
                                   p=p,
                                   diff_stencil_arm=diff_stencil_arm)
cp2, distance2, _, _ = s.closest_point(grid)
assert np.allclose(cp, cp2)
assert np.allclose(distance, distance2)

# Corners of the virtual grid, superset of `grid`
ll = np.array(3 * [grid.min()]) - 3 * dx
ur = np.array(3 * [grid.max()]) + 3 * dx
virtual_grid_shape = np.abs(ur-ll) / dx + 1

# The (i,j,...) indices of the grid points, taking `ll` as origin.
int_grid = np.round((grid - ll) / dx).astype(np.int)

# Initial conditions
th, phi, r = cart2sph(grid[:, 0], grid[:, 1], grid[:, 2])
u = np.cos(phi + np.pi / 2)
# Let's keep a copy of the initial conditions
initial_u = u.copy()

# Build interpolation and differential matrix.
E = build_interp_matrix(int_grid, cp, dx, p, ll, virtual_grid_shape)
L = build_diff_matrix(int_grid, dx, virtual_grid_shape)


xp, yp, zp = s.parametric_grid(65)
_, phi_plot, _ = cart2sph(xp, yp, zp)
Eplot = build_interp_matrix(int_grid,
                            np.column_stack((xp.ravel(),
                                             yp.ravel(),
                                             zp.ravel())),
                            dx, p, ll, virtual_grid_shape)

if PLOT:
    # Plotting code. Build a pipeline to be able to change the data later.
    src = mlab.pipeline.grid_source(xp, yp, zp,
                                    scalars=(Eplot * u).reshape(xp.shape))
    normals = mlab.pipeline.poly_data_normals(src)
    surf = mlab.pipeline.surface(normals)
    mlab.colorbar()


Tf = 2
dt = 0.1 * np.min(dx)**2
numtimesteps = int(Tf // dt + 1)
errors = []  # To store the error at each timestep
# Explicit Forward Euler time stepping
for kt in xrange(numtimesteps):
    unew = u + dt * (L*u)
    u = E*unew
    t = kt * dt
    if not kt%100 or kt == (numtimesteps-1):
        print "time: {0:2f}, {1:2f} %".format(t, float(kt) / numtimesteps)
        sphplot = Eplot * u
        true_solution = np.exp(-2*t) * np.cos(phi_plot + np.pi / 2)
        step_error = (np.abs(true_solution - sphplot.reshape(xp.shape)).sum() /
                      np.abs(true_solution).sum())
        errors.append(step_error)
        if PLOT:
            src.data.point_data.scalars = sphplot
            src.data.point_data.scalars.name = 'scalars'
            src.data.modified()