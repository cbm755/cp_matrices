"""Solves the PDE problems on the surface of a brain.
Does not output PETSc matrices.
"""
import numpy as np
import scipy as sp
import pickle
import math


from cp.surfaces import Mesh
from cp.tools.io import load_ply
from cp.build_matrices import build_interp_matrix, build_diff_matrix

PLOT = False
OUTPUT_PETSc = True

if PLOT:
    try:
        from mayavi import mlab
    except ImportError:
        from enthought.mayavi import mlab

# output options
basename = 'brain_r777'

# Load vertices and faces, and instantiate surface
plyscale = 0;
vert, faces = load_ply('brain-lh_scale_' + str(plyscale) + '.ply')
m = Mesh(vert, faces)

p = 3
diff_stencil_arm = 1
dim = 3

print 'building grid'
index, distance, grid, dx = m.grid(num_blocks_per_dim=101,
                                   levels=1,
                                   p=p,
                                   diff_stencil_arm=diff_stencil_arm)
print 'finding closest points'
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
u = np.zeros(grid[:,0].size)
# Let's keep a copy of the initial conditions
initial_u = u.copy()

# Let's build the matrices. TODO: I think it would be nicer to use
# `grid` instead of `int_grid`. It is a simple change.
print 'building E matrix'
E = build_interp_matrix(int_grid, cp, dx, p, ll, virtual_grid_shape)
# TODO: being able to select different laplacian operators. Currently
# it uses the second order laplacian for 2D and 3D. We could probably
# use stencils.py, and give the stencil as a parameter
print 'building L matrix'
L = build_diff_matrix(int_grid, dx, virtual_grid_shape)

# Points in the surface of the sphere, used por plotting
#xp, yp, zp = Sphere().parametric_grid(65)
#_, phi_plot, _ = cart2sph(xp, yp, zp)
print 'building Eplot matrix'
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
    surf.module_manager.scalar_lut_manager.use_default_range = False
    surf.module_manager.scalar_lut_manager.data_range = (0.25, 0.45)


# reaction-diffusion equation
# parameters:
alpha = 50.0     # (model) coefficient of reaction term
gammaS = 100.0      # (numerical) coefficient of sources term
v0 = 0.5          # magnitude of point sources

# load source locations
sources = np.loadtxt("brain_sources.txt")
#sources = pickle.load(file('brain_sources.pickle'))
nsrcs = sources.shape[0]    # number of sources
# sourcecp, _, _, _ = m.closest_point(sources)
# build the source term
v = 0
varsq = 5*dx[0]      # scale delta fns somehow
# second run: (delete a point)
#for srccount in xrange(nsrcs-1):
for srccount in xrange(nsrcs):
  vdist = (grid[:,0]-sources[srccount,0])**2 + (grid[:,1]-sources[srccount,1])**2 + (grid[:,2]-sources[srccount,2])**2
  # exp fns around sources
  v = v + 1.0/math.sqrt(2*math.pi*varsq) * np.exp( -vdist / (2 * varsq))
# cp-ext
v = E*v

turnoff = 1
if PLOT:
    mlab.points3d(sources[:,0], sources[:,1], sources[:,2], scale_factor=0.1,color=(1,1,1))
    #mlab.points3d(sources[turnoff,0], sources[turnoff,1], sources[turnoff,2], scale_factor=0.1, color=(1,1,1))
    #raw_input("press enter to continue")

# choose a Closest Point Method algorithm
cpm = 1

# choose timestep
if cpm == 0:
    dt = 0.2 * np.min(dx)**2
elif cpm == 1:
    dt = 0.2 * np.min(dx)**2
elif cpm == 2:
    dt = 0.5 * np.min(dx)

Tf = 0.5
numtimesteps = int(np.ceil(Tf / dt))
turn_off_at_time = 10.0
turn_off_at = int(np.ceil(turn_off_at_time / dt))
dt = Tf / numtimesteps
print "will turn off some sources at kt=" + str(turn_off_at)

# build the vGMM matrix
print 'assemble vGMM matrix'
if cpm == 1 or cpm == 2:
    #I = speye(L.shape[0], L.shape[1])
    I = sp.sparse.eye(*L.shape)
    lamb = 4.0/np.min(dx)**2
    M = E*L - lamb*(I - E)
if cpm == 2:
    A = I - dt*M


print 'starting time-stepping'
# Explicit Forward Euler time stepping
for kt in xrange(numtimesteps):
    if kt == turn_off_at:
        # TODO: should  be a funciton or something
        print 'turning some srcs off'
        v = 0
        for srccount in xrange(nsrcs):
            # check whether this source is in the list of those to be turned off
            if srccount in (turnoff,):
                if PLOT:
                    mlab.points3d(sources[srccount,0], sources[srccount,1], sources[srccount,2], scale_factor=0.1,color=(0,0,0))
            else:
                vdist = (grid[:,0]-sources[srccount,0])**2 + (grid[:,1]-sources[srccount,1])**2 + (grid[:,2]-sources[srccount,2])**2
                # exp fns around sources
                v = v + np.exp( -vdist / (2 * varsq))
        # cp-ext
        v = E*v


    if cpm == 0:
        # explicit Euler, Ruuth--Merriman style
        unew = u + dt * (L*u)
        u = E*unew
    elif cpm == 1:
        # explicit Euler, von Glehn--Maerz--Macdonald

        # TODO: clean this up!
        # this now solves a nonlinear reaction-diffusion equation
        unew = u + dt * (M*u - alpha*u/(1+u) - gammaS*(u*v - v0*v) )
        u = unew;
    elif cpm == 2:
        # implicit Euler, vGMM
        unew = splinalg.spsolve(A, u)
        u = unew;
    # unew = u + dt * (L*u)
    # u = E*unew
    t = kt * dt
    if not kt%25 or kt == (numtimesteps-1):


        print "kt={:d}, time: {:2f}, {:2f} %, u_minmax = [{:g},{:g}]".format(kt, t, 100 * float(kt) / numtimesteps, u.min(), u.max())
        uplot = Eplot * u

        # output data
        fname = '{:s}_gridsoln_kt_{:0>6d}.pickle'.format(basename, kt)
        pickle.dump((u), file(fname, 'w'))
        # more output, as binary float32 data:
        fname = '{:s}_plot_scale{:d}_kt_{:0>6d}.ascii'.format(basename, plyscale, kt)
        uplot.astype('f').tofile(fname)
        # or ascii output
        uplot.astype('f').tofile(fname, sep=" ", format="%s")

        #true_solution = np.exp(-2*t) * np.cos(phi_plot + np.pi / 2)
        #step_error = (np.abs(true_solution - sphplot.reshape(xp.shape)).sum() /
        #np.abs(true_solution).sum())
        #errors.append(step_error)
        if PLOT and not kt%100:
            src.data.point_data.scalars = uplot
            src.data.point_data.scalars.name = 'scalars'
            src.data.modified()
            raw_input("press enter to continue")


if OUTPUT_PETSc:
    print 'saving matrices to petsc format on disk'
    import cp.tools.scipy_petsc_conversions as conv
    st = timeit.default_timer()
    if cpm == 0:
        conv.save_scipy_to_petsc_ondisk(L, (7,0), 'brain_Lmatrix.dat')
        conv.save_scipy_to_petsc_ondisk(E, (64,0), 'brain_Ematrix.dat')
    elif cpm == 1:
        # 2.5*64 + small safety factor
        conv.save_scipy_to_petsc_ondisk(M, (165,0), 'brain_Mmatrix.dat')
    elif cpm == 2:
        conv.save_scipy_to_petsc_ondisk(A, (165,0), 'brain_Amatrix.dat')

    final_u = u
    print 'saving dx, ICs, soln to disk'
    # todo: v not right yet
    pickle.dump((dx, cpm, Tf, numtimesteps, dt, initial_u, final_u, alpha, gammaS, v0, v), file('brain_nonpetsc_data.pickle','w'))
    timeit.default_timer() - st
