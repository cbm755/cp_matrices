"""Solves the heat equation on a true sphere."""
import numpy as np
from scipy.sparse import eye
import matplotlib.pyplot as plt
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab

from cp.surfaces import Sphere
from cp.build_matrices import build_interp_matrix, build_diff_matrix
from cp.surfaces.coordinate_transform import cart2sph, cart2pol


PLOT = True

s = Sphere()

p = 3
diff_stencil_arm = 1
dim = 3

# As a byproduct of finding the banded grid, we already have its
# closest points, so we don't really have to call s.closest_point()
cp, distance, grid, dx = s.grid(num_blocks_per_dim=41,
                                   levels=2,
                                   p=p,
                                   diff_stencil_arm=diff_stencil_arm)

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
I = eye(*E.shape)

# forcing factor
gslam = -0.05 / 5
# plane of image projection
surf_plane = 0

# Load png and select non-zero channel
image = plt.imread('examples/RD_mask.png')[..., 3]
vertical_res, horizontal_res = image.shape

# Image limits
xmin, xmax = -1.5, 1.5
ymin, ymax = -1.5, 1.5
xi = np.round(horizontal_res * (grid[:, 0] - xmin) / (xmax - xmin)).astype(np.int)
yi = (vertical_res - np.round(vertical_res * (grid[:, 1] - ymin) /
                             (ymax - ymin))).astype(np.int)
# Forcing functions
chi = np.zeros(grid.shape[0])
# Create forcing function
ix = grid[:, 2] >= surf_plane
chi[ix] = image[xi[ix], yi[ix]]
# make sure forcing is a cp extension
chi = E * chi
# re thresold
chi = (chi > 0.9).astype(np.int)

xp, yp, zp = s.parametric_grid(256)
_, phi_plot, _ = cart2sph(xp, yp, zp)
Eplot = build_interp_matrix(int_grid,
                            np.column_stack((xp.ravel(),
                                             yp.ravel(),
                                             zp.ravel())),
                            dx, p, ll, virtual_grid_shape)

# Plot forcing function
mlab.figure(1)
mlab.mesh(xp, yp, zp, scalars=(Eplot*(1-chi)).reshape(xp.shape))
mlab.title('forcing function')

# parameters and functions for Gray--Scott
# 120 works with 0.025
F = 0.054
kk = 0.063
nuu = 1 / (3/dx.min())**2
nuv = nuu / 2.
f = lambda u, v: -u * v**2 + F * (1 - u)
g = lambda u, v: u * v**2 - (F+kk) * v

# Initial conditions - small perturbation from steady state
th, r, _ = cart2pol(*grid.T)
pert = 0.5 * np.exp(-(10*(grid[:,2] - 1))**2) + 0.5 * np.random.randn(grid.shape[0])
u0 = 1 - pert
v0 = 0 + 0.5 * pert
u = u0
v = v0

if PLOT:
    mlab.figure(2)
    # Plotting code. Build a pipeline to be able to change the data later.
    src = mlab.pipeline.grid_source(xp, yp, zp,
                                    scalars=np.zeros(xp.shape))
    normals = mlab.pipeline.poly_data_normals(src)
    surf = mlab.pipeline.surface(normals)
    mlab.colorbar()


# cpmol
lambda_ = 4 * nuu / (dx.min()**2)
Au = nuu * E * L - lambda_ * (I - E)
Av = nuv * E * L - lambda_ * (I - E)

Tf = 5000
dt = 0.1 * (1. / max(nuu, nuv)) * np.min(dx)**2
print "dt", dt
numtimesteps = int(Tf // dt + 1)
# Explicit Forward Euler time stepping
for kt in xrange(numtimesteps):
    if kt < 10:  # No forcing
        unew = u + dt * (E * f(u,v) + Au * u)
        vnew = v + dt * (E * g(u, v) + Av * v)
    elif kt < 3000:  # Turn on first forcing function
        unew = u + dt * (E * f(u, v) + Au*u + gslam * chi * (u-0.3))
        vnew = v + dt * (E * g(u, v) + Av*v + gslam * chi * (v-0.6))
    u = unew.copy()
    v = vnew.copy()
    t = kt * dt
    if not kt%20 or kt == (numtimesteps-1) or kt <=50:
        print "time: {0:2f}, {1:2f} %".format(t, float(kt) / numtimesteps)
        sphplot = Eplot * u
        if PLOT:
            src.data.point_data.scalars = sphplot
            src.data.point_data.scalars.name = 'scalars'
            src.data.modified()