import numpy as np

from cp.surfaces import coordinate_transform
from cp.tools.io import load_ply
from cp.mesh_surface import Mesh
from cp.build_matrices import build_interp_matrix, build_diff_matrix
from cp.surfaces.coordinate_transform import cart2sph, sph2cart


v, f = load_ply('cp/tests/data/sphere_refined.ply')
m = Mesh(v, f)
# If I set levels=3, somehow points in the evolution band (ie,
# basepoints+interpolation_offsets) end up outside the band :S Doesn't
# make any sense...
index, grid, distance, dx = m.grid(num_blocks=41, levels=1)
cp, dist, _, _ = m.closest_point(index, grid)
th, phi, r = cart2sph(*grid.T)
cpx, cpy, cpz = sph2cart(th, phi, 1.)
cp_real = np.column_stack((cpx, cpy, cpz))
dist_real = np.abs(r - 1.)

ll = np.array(3 * [grid.min()]) - 3 * dx
ur = np.array(3 * [grid.max()]) + 3 * dx
virtual_grid_shape = np.abs(ur-ll) / dx + 1

int_grid = np.round((grid - ll) / dx).astype(np.int)

dim = 3
p = 3
order = 2

th, phi, r = coordinate_transform.cart2sph(grid[:, 0], grid[:, 1], grid[:, 2])
u = np.cos(phi + np.pi / 2)

initial_u = u.copy()

E = build_interp_matrix(int_grid, cp, dx, p, ll, virtual_grid_shape)
L = build_diff_matrix(int_grid, dx, virtual_grid_shape)

###
Tf = 1
dt = 0.1 * np.min(dx)**2
numtimesteps = int(Tf // dt + 1)
err = []
#th_plot, phi_plot, r = cart2sph(*cp.T)
th_plot, phi_plot = np.mgrid[-np.pi:np.pi:65j,
                             -np.pi/2:np.pi/2:65j].reshape((2, -1))
xp, yp, zp = sph2cart(th_plot, phi_plot, 1)
Eplot = build_interp_matrix(int_grid, np.column_stack((xp, yp, zp)), dx, p, ll, virtual_grid_shape)
from mayavi import mlab
s = mlab.points3d(xp, yp, zp, Eplot * u)
mlab.colorbar()

for kt in xrange(numtimesteps):
    unew = u + dt * (L*u)
    u = E*unew
    t = kt * dt
    if not kt%100:
        print round(float(kt) / numtimesteps, 2)
        sphplot = Eplot * u
        true_solution = np.exp(-2*t) * np.cos(phi_plot + np.pi / 2)
        err.append(np.abs(true_solution - sphplot).sum() / np.abs(true_solution).sum())
        s.mlab_source.scalars = sphplot
        
