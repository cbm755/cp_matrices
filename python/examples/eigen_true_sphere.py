"""Solves the heat equation on a true sphere."""
from time import time
import numpy as np
import scipy.sparse.linalg as splinalg
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab

from cp.surfaces import Sphere
from cp.build_matrices import build_interp_matrix, build_diff_matrix, build_linear_diagonal_splitting
from cp.surfaces.coordinate_transform import cart2sph


PLOT = True

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
ll = np.array(dim * [grid.min()]) - 3 * dx
ur = np.array(dim * [grid.max()]) + 3 * dx
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
M = build_linear_diagonal_splitting(L, E)

# Compute eigenvalues
t = time()
print "Computing eigenvalues..."
Evals, Evecs = splinalg.eigs(-M, k=32, which="SM")
sorted_indices = np.argsort(Evals)
print "...took", time() -t
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
                                    scalars=zp)
    normals = mlab.pipeline.poly_data_normals(src)
    surf = mlab.pipeline.surface(normals)
    mlab.colorbar()


for ii in xrange(9, 15+1):
    evec = Evecs[:, sorted_indices[ii]]
    emax = evec.real.max()
    emin = evec.real.min()
    absmax = max(-emin, emax)
    src.data.point_data.scalars = Eplot * E * evec.real
    src.data.point_data.scalars.name = 'scalars'
    src.data.modified()
