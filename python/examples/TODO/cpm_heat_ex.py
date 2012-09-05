r"""
Solve the diffusion (heat) equation

     $u_t = \kappa \triangle_s u$

on a surface.  $\triangle_s$ is the Laplace--Beltrami operator.
$\kappa$ is a scalar coefficient.

Uses explicit forward Euler timestepping.
"""
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import  mlab
from time import time
from operator import mod


usz = M.shape[0]  # system size
# Tfinal = 1
Nsteps = 500      # How many time steps
dt = 0.1*dx*dx    # note O(dx^2) b/c this is explicit
kappa = 1.0       # diffusion coefficient


# setup initial conditions
u0 = 1*np.random.randn(usz)
# scale to max abs 1
u0 = u0 / max(abs(u0.max()), abs(u0.min()))
u0plot = Eplot*E*u0
u = u0


# setup viz
f = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(640,640))
mlab.clf()
if (1==1):
    # build a pipeline so we can explicitly change the data later
    src = mlab.pipeline.grid_source(x, y, z, scalars=u0plot.reshape(x.shape))
    normals = mlab.pipeline.poly_data_normals(src)
    surf = mlab.pipeline.surface(normals)
else:
    # this is easier but less convenient to change the data
    s = mlab.mesh(x, y, z, scalars=real(uplot.reshape(x.shape)))
# TODO: size is not fontsize
mlab.title('Initial conditions', size=0.2)
mlab.show()
pause()


# Explicit Forward Euler time stepping
for kt in range(0,Nsteps):
    now = time()
    utilde = E*u
    unew = utilde[0:usz] + (dt*kappa)*(D*utilde)
    print "timestep " + str(kt) + ", took %.3g s" % (time()-now)

    t = 0 + dt*(kt+1)

    # Do viz for some timesteps
    if (kt < 10) or (mod(kt+1,20) == 0):
        uplot = Eplot*(E*unew)  # brackets necessary for speed here
        # update the viz
        src.data.point_data.scalars = Eplot*(E*unew)
        src.data.point_data.scalars.name = 'scalars'
        src.data.modified()
        mlab.title("time = " + str(t) + ", step #" + str(kt+1), size=0.2)
        mlab.show()
        #mlab.savefig('hemisphere' + str(ii) + '_' + str(eval) + '.png')

        # A simpler but less efficient approach: draw a new surface each time
        if (1==0):
            mlab.clf()
            s = mlab.mesh(x, y, z, scalars=uplot.reshape(x.shape))
            #, vmin=-absmax, vmax=absmax)
            mlab.title("time = " + str(t) + ", step #" + str(kt+1))
            mlab.show()
            pause(0.1)

    u = unew

