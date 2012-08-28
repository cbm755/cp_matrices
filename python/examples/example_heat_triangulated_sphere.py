import numpy as np

from cp.surfaces import coordinate_transform
from cp.tools.io import load_ply
from cp.mesh_surface import Mesh
from cp.build_matrices import build_interp_matrix, build_diff_matrix


v, f = load_ply('cp/tests/data/sphere_refined.ply')
m = Mesh(v, f)
# If I set levels=3, somehow points in the evolution band (ie,
# basepoints+interpolation_offsets) end up outside the band :S Doesn't
# make any sense...
index, grid, distance, dx = m.grid(num_blocks=41, levels=1)
cp, dist, _, _ = m.closest_point(index, grid)

ll = np.array(3 * [grid.min()]) - 3 * dx
ur = np.array(3 * [grid.max()]) + 3 * dx
virtual_grid_shape = np.abs(ur-ll) / dx + 1

int_grid = np.round((grid - ll) / dx).astype(np.int)

dim = 3
p = 3
order = 2

th, phi, r = coordinate_transform.cart2sph(grid[:, 0], grid[:, 1], grid[:, 2])
u = np.cos(phi + np.pi / 2) + 1

initial_u = u.copy()

E = build_interp_matrix(int_grid, cp, dx, p, ll, virtual_grid_shape)
L = build_diff_matrix(int_grid, dx, virtual_grid_shape)

###
from mayavi import mlab
s = mlab.points3d(cp[:,0], cp[:,1], cp[:,2], initial_u)
mlab.colorbar()
Tf = 1
dt = 0.1 * dx[0]**2
numtimesteps = int(Tf // dt + 1)
l = []
for kt in xrange(numtimesteps):
    unew = u + dt * (L*u)
    u = E*unew
    t = kt * dt
    l.append(u.mean())
    if not kt%100:
        s.mlab_source.scalars = u
