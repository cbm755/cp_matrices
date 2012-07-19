r"""
Solve the brusselator model for pattern formation on a surface.

The system is two coupled PDEs for concentrations u and v.  These
diffuse linearly and interact with each other nonlinearly.

     $u_t = f(u,v) + \nu_u \triangle_s u$
     $v_t = g(u,v) + \nu_v \triangle_s v$

$\triangle_s$ is the Laplace--Beltrami operator.

$\nu_u$ and $\nu_v$ are scalar diffusion coefficients.

Uses implicit ebdf2 timestepping.
"""

from scipy.sparse import identity
#from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import gmres as gmres
#from scipy.linsolve import spsolve
#import scipy.sparse.linalg.dsolve as spsolve
from scipy.sparse.linalg import spsolve
import numpy as np

try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab
from time import time
from operators import mod


# Parameters
a = 3.0
b = 10.2

## stripy patterns
scale = 30**2  # controls the spatial size of the patterns
nuu = 3.8/(scale)
## spotty patterns
#scale = 32**2
#nuu = 2.5/(scale)
## holey patterns
#scale = 26**2
#nuu = 5.0/(scale)

nuv = 10.0/(scale)

def f(u,v):
    return a - (b+1)*u + u*u*v
def g(u,v):
    return b*u - u*u*v
gmtol = 1e-9





usz = M.shape[0]  # system size
# Tfinal = 1
Nsteps = 1000      # How many time steps
dt = 0.5*dx       # note O(dx) b/c this is implicit
Iterative = True  # using  an iterative solver is faster


I = identity(usz)

# Setup initial conditions (random perturbations around the
# zero-diffusion steady state)
# TODO: look at the scaling 0.1 here...
u0 = a*np.ones(usz) + 0.1*np.random.randn(usz)
v0 = b/a*np.ones(usz) + 0.1*np.random.randn(usz)
if (1==0):  # smooth with one step of heat equation
    A = I - (0.1*dx)*M
    (u0,flag1) = gmres(A, u0, x0=u0, tol=1e-5, restrt=10)
    (v0,flag2) = gmres(A, v0, x0=v0, tol=1e-5, restrt=10)
    A = None
u0plot = Eplot*E*u0
v0plot = Eplot*E*v0
u = u0
v = v0



# setup viz
fig = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(640,640))
mlab.clf()
# build a pipeline so we can explicitly change the data later
src = mlab.pipeline.grid_source(x, y, z, scalars=u0plot.reshape(x.shape))
normals = mlab.pipeline.poly_data_normals(src)
surf = mlab.pipeline.surface(normals)
# TODO: size is not fontsize
mlab.title('u: initial conditions', size=0.2)
mlab.show()
pause(1)



# Build the backward Euler time-stepping matrices
Au_be = I - (nuu*dt)*M
Av_be = I - (nuv*dt)*M

# timestepping
for kt in range(0,Nsteps):
    now = time()
    if (kt == 0):
        (unew,flag1) = gmres(Au_be, u + dt*f(u,v), x0=u, tol=gmtol, restrt=10)
        (vnew,flag2) = gmres(Av_be, v + dt*g(u,v), x0=v, tol=gmtol, restrt=10)
        # don't need the BE matrices anymore
        Au_be = None;  Av_be = None
        # build new bdf2 matrices for the remaining timesteps
        Au_bdf2 = I - (nuu*2.0*dt/3.0)*M
        Av_bdf2 = I - (nuv*2.0*dt/3.0)*M
    else:
        (unew,flag1) = gmres(Au_bdf2, \
                                 4.0/3.0*u - 1.0/3.0*uold + \
                                 4.0*dt/3.0*f(u,v) - 2.0*dt/3.0*f(uold,vold), \
                                 x0=u, tol=gmtol)
        (vnew,flag2) = gmres(Av_bdf2, \
                                 4.0/3.0*v - 1.0/3.0*vold + \
                                 4.0*dt/3.0*g(u,v) - 2.0*dt/3.0*g(uold,vold), \
                                 x0=v, tol=gmtol)
    if ((flag1 != 0) or (flag2 != 0)):
        print "  flags = " + str((flag1,flag2))

    print "timestep " + str(kt) + ", took %.3g s" % (time()-now)

    t = 0 + dt*(kt+1)

    # Do viz for some timesteps
    if (kt < 10) or (mod(kt+1,25) == 0):
        uplot = Eplot*(E*unew)  # brackets necessary for speed here
        # update the viz
        src.data.point_data.scalars = Eplot*(E*unew)
        src.data.point_data.scalars.name = 'scalars'
        src.data.modified()
        mlab.title("time = " + str(t) + ", step #" + str(kt+1), size=0.2)
        mlab.show()
        #mlab.savefig('hemisphere' + str(ii) + '_' + str(eval) + '.png')

    uold = u
    vold = v
    u = unew
    v = vnew

