"""Solves a Fitzhugh-Nagumo problem on the surface of a brain.
Outputs binary files of the solution vector
"""
import numpy as np
import scipy as sp
import pickle
import math


from cp.surfaces import Mesh
from cp.tools.io import load_ply
from cp.build_matrices import build_interp_matrix, build_diff_matrix
#from cp.surfaces.coordinate_transform import cart2sph

PLOT = True

if PLOT:
    try:
        from mayavi import mlab
    except ImportError:
        from enthought.mayavi import mlab

# output options
basename = 'brain_fhn001'

# Load vertices and faces, and instantiate surface
plyscale = 0;
vert, faces = load_ply('brain-lh_scale_' + str(plyscale) + '.ply')
m = Mesh(vert, faces)

p = 3
diff_stencil_arm = 1
dim = 3

index, distance, grid, dx = m.grid(num_blocks_per_dim=160,
                                   levels=1,
                                   p=p,
                                   diff_stencil_arm=diff_stencil_arm)
cp, dist, _, _ = m.closest_point(grid, index)

# The points in `grid` can be thought to be a subset of a virtual grid
# (for instance, the result of meshgrid). `ll` is the lower left
# corner of such virtual grid, and `ur` the upper right corner. The
# padding (\pm 3 * dx) should not be needed, but I haven't check if
# that's true. Fixing the shape of this virtual grid let's us easily
# go from an (i,j,...) index to a linear index and back.
ll = np.array(dim * [grid.min()]) - 3 * dx
ur = np.array(dim * [grid.max()]) + 3 * dx
virtual_grid_shape = np.abs(ur-ll) / dx + 1

# The (i,j,...) indices of the grid points, taking `ll` as origin.
int_grid = np.round((grid - ll) / dx).astype(np.int)

# To set the initial conditions I directly set a value for each grid
# point. Another option would be to set the values in the vertices,
# then interpolate to get the values in each closest point (using
# scipy.interpolate.griddata) and finally extend that to the whole
# grid
#th, phi, r = cart2sph(grid[:, 0], grid[:, 1], grid[:, 2])
#u = np.cos(phi + np.pi / 2)
#u = grid[:,0]**3
u = (grid[:,0]>0)&(grid[:,1]>0)&(grid[:,2]>0);
v = (grid[:,0]<0)&(grid[:,1]>0)&(grid[:,2]>0);
# Let's keep a copy of the initial conditions
initial_u = u.copy()
initial_v = v.copy()

# Let's build the matrices. TODO: I think it would be nicer to use
# `grid` instead of `int_grid`. It is a simple change.
E = build_interp_matrix(int_grid, cp, dx, p, ll, virtual_grid_shape)
# TODO: being able to select different laplacian operators. Currently
# it uses the second order laplacian for 2D and 3D. We could probably
# use stencils.py, and give the stencil as a parameter
L = build_diff_matrix(int_grid, dx, virtual_grid_shape)

# Points in the surface of the sphere, used por plotting
#xp, yp, zp = Sphere().parametric_grid(65)
#_, phi_plot, _ = cart2sph(xp, yp, zp)
Eplot = build_interp_matrix(int_grid,
                            m.vertices,
                            dx, 1, ll, virtual_grid_shape)

if PLOT:
    # Plotting code. Build a pipeline to be able to change the data later.
    #src = mlab.pipeline.grid_source(xp, yp, zp,
    #                                scalars=(Eplot * u).reshape(xp.shape))
    src = mlab.pipeline.triangular_mesh_source(m.vertices[:,0],
                                               m.vertices[:,1],
                                               m.vertices[:,2], m.faces,
                                               scalars=(Eplot * u))
    normals = mlab.pipeline.poly_data_normals(src)
    surf = mlab.pipeline.surface(normals)
    mlab.colorbar()


# reaction-diffusion equation (Fitzhugh-Nagumo)
# parameters:
a = 0.1
eps = 0.015
beta = 0.5
gamma = 1.0
delta = 0

nuu = 1.44e-4
nuv = 3.6e-6
f = lambda u, v: (a-u) * (u-1) * u - v
g = lambda u, v: eps * (beta*u - gamma*v - delta)

# choose a Closest Point Method algorithm
cpm = 1

# choose timestep
if cpm == 0:
    dt = 0.2 * (1. / max(nuu, nuv)) * np.min(dx)**2
elif cpm == 1:
    dt = 0.2 * (1. / max(nuu, nuv)) * np.min(dx)**2
elif cpm == 2:
    dt = 0.5 * (1. / max(nuu, nuv)) * np.min(dx)

Tf = 2000.
numtimesteps = int(np.ceil(Tf / dt))
dt = Tf / numtimesteps


# build the vGMM matrix
if cpm == 1 or cpm == 2:
    #I = speye(L.shape[0], L.shape[1])
    I = sp.sparse.eye(*L.shape)
    lamb = max(nuu,nuv) * 4.0/np.min(dx)**2   # TODO: should this also be scaled to nuu?
    Mu = nuu * E*L - lamb*(I - E)
    Mv = nuv * E*L - lamb*(I - E)

if cpm == 2:
    print "Warning: implicit not yet implemented"
    A = I - dt*Mu - dt*Mv



errors = []  # To store the error at each timestep
# Explicit Forward Euler time stepping
for kt in xrange(numtimesteps):

    if cpm == 0:
        # explicit Euler, Ruuth--Merriman style
        unew = u + dt * (f(u,v) + nuu * L * u)
        vnew = v + dt * (g(u,v) + nuv * L * v)
        u = E * unew
        v = E * vnew
    elif cpm == 1:
        # explicit Euler, von Glehn--Maerz--Macdonald
        unew = u + dt * (E * f(u,v) + Mu * u)
        vnew = v + dt * (E * g(u,v) + Mv * v)
        u = unew;
        v = vnew;
    elif cpm == 2:
        # implicit Euler, vGMM
        unew = splinalg.spsolve(A, u)
        u = unew;

    t = kt * dt
    if not kt%5 or kt == (numtimesteps-1):


        print "kt={:d}, time: {:2f}, {:2f} %, u_minmax = [{:g},{:g}]".format(kt, t, 100 * float(kt) / numtimesteps, u.min(), u.max())
        uplot = Eplot * u

        # output data
        fname = '{:s}_gridsoln_kt_{:0>6d}.pickle'.format(basename, kt)
        pickle.dump((u), file(fname, 'w'))
        # more output, as binary float32 data:
        fname = '{:s}_plot_scale{:d}_kt_{:0>6d}.bin'.format(basename, plyscale, kt)
        uplot.astype('f').tofile(fname)

        #true_solution = np.exp(-2*t) * np.cos(phi_plot + np.pi / 2)
        #step_error = (np.abs(true_solution - sphplot.reshape(xp.shape)).sum() /
        #np.abs(true_solution).sum())
        #errors.append(step_error)
        if PLOT:
            src.data.point_data.scalars = uplot
            src.data.point_data.scalars.name = 'scalars'
            src.data.modified()
            raw_input("press enter to continue")

