"""Solves the heat equation on a true sphere."""
import numpy as np
import pickle
import timeit

from scipy.sparse import eye as speye
import scipy.sparse.linalg as splinalg

from cp.surfaces import Circle
from cp.build_matrices import build_interp_matrix, build_diff_matrix
from cp.surfaces.coordinate_transform import cart2pol

s = Circle()

p = 3
diff_stencil_arm = 1
dim = 2

# As a byproduct of finding the banded grid, we already have its
# closest points, so we don't really have to call s.closest_point()
cp, distance, grid, dx = s.grid(num_blocks_per_dim=321,
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
th, r = cart2pol(grid[:, 0], grid[:, 1])
u = np.cos(th + np.pi / 2)
# Let's keep a copy of the initial conditions
initial_u = u.copy()

# Build interpolation and differential matrix.
E = build_interp_matrix(int_grid, cp, dx, p, ll, virtual_grid_shape)
L = build_diff_matrix(int_grid, dx, virtual_grid_shape)


xp, yp = s.parametric_grid(65)
th_plot, _ = cart2pol(xp, yp)
Eplot = build_interp_matrix(int_grid,
                            np.column_stack((xp.ravel(),
                                             yp.ravel())),
                            dx, p, ll, virtual_grid_shape)


Tf = 2
use_implicit = True
if use_implicit:
    dt = 0.5 * np.min(dx)
    #I = speye(L.shape[0], L.shape[1])
    I = speye(*L.shape)
    lamb = 4.0/np.min(dx)**2
    M = E*L - lamb*(I - E)
    A = I - dt*M
else:
    dt = 0.1 * np.min(dx)**2

numtimesteps = int(Tf // dt + 1)

start_time = timeit.default_timer()

# Explicit Forward Euler time stepping
for kt in xrange(numtimesteps):
    # explicit euler
    #unew = u + dt * (L*u)
    #u = E*unew
    # explicit MOL
    #unew = u + dt * (M*u)
    # implicit MOL
    unew = splinalg.spsolve(A, u)

    t = kt * dt
    if not kt%10 or kt == (numtimesteps-1):
        uplot = Eplot * u
        true_solution = np.exp(-t) * np.cos(th_plot + np.pi / 2)
        max_error = (np.abs(true_solution - uplot)).max()
        print "time: {0:2f} ({1:2.2f}%), err={2:g}".format(t, 100*float(kt) / numtimesteps, max_error)
        #errors.append(step_error)
        #if PLOT:
        #    src.data.point_data.scalars = sphplot
        #    src.data.point_data.scalars.name = 'scalars'
        #    src.data.modified()

print "Serial code time=", timeit.default_timer() - start_time


print 'saving matrices to petsc format on disk'
import cp.tools.scipy_petsc_conversions as conv
conv.save_scipy_to_petsc_ondisk(L, 'Lmatrix.dat')
conv.save_scipy_to_petsc_ondisk(E, 'Ematrix.dat')
conv.save_scipy_to_petsc_ondisk(A, 'Amatrix.dat')

final_u = u
print 'saving dx, ICs, soln to disk'
pickle.dump((dx,initial_u, final_u), file('non_petsc_data.pickle','w'))
