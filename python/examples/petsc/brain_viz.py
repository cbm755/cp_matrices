
"""Solves the heat equation on a triangulated sphere."""
import numpy as np
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab

import pickle
from cp.surfaces import Mesh
# Since our mesh is a sphere, we'll take advantage of its
# parametric_plot method
from cp.tools.io import load_ply
#from cp.build_matrices import build_interp_matrix, build_diff_matrix
# TODO: move coordinate_transform out of cp.surfaces (maybe to
# cp.tools?)
#from cp.surfaces.coordinate_transform import cart2sph



# which input to read?
basename = 'brain_r001'

# Load vertices and faces, and instantiate surface
plyscale = 0;
vert, faces = load_ply('brain-lh_scale_' + str(plyscale) + '.ply')
#m = Mesh(vert, faces)


if (1==1):
    # Plotting code. Build a pipeline to be able to change the data later.
    #src = mlab.pipeline.grid_source(xp, yp, zp,
    #                                scalars=(Eplot * u).reshape(xp.shape))
    src = mlab.pipeline.triangular_mesh_source(vert[:,0], vert[:,1], vert[:,2],
            faces, scalars=np.ones((vert.shape[0],)))

    normals = mlab.pipeline.poly_data_normals(src)
    surf = mlab.pipeline.surface(normals)
    mlab.colorbar()

#Tf = 0.2
#dt = 0.2 * np.min(dx)**2
#numtimesteps = int(Tf // dt + 1)
#dt = Tf / numtimesteps

mlab.options.offscreen = True

if (1==0):
    #(dx, Tf, numtimesteps, dt, initial_u_plot) = pickle.load(file('brain_IC_plotvec.pickle'))
    #(dx, Tf, numtimesteps, dt, uplot) = pickle.load(file('brain_final_plotvec.pickle'))
    uplot = np.fromfile('testout.bin', 'f')
    src.data.point_data.scalars = uplot
    src.data.point_data.scalars.name = 'scalars'
    src.data.modified()
    raw_input("press enter to continue")


# or load from the binary files:
for kt in xrange(0, 3100, 300):
    fname = '{:s}_plot_scale{:d}_kt_{:0>6d}.bin'.format(basename, plyscale, kt)
    uplot = np.fromfile(fname, 'f')

    # TODO: fix lower/upper bounds
    src.data.point_data.scalars = uplot
    src.data.point_data.scalars.name = 'scalars'
    src.data.modified()
    # TODO: how to do pause(0) or drawnow() in mayavi?
    print 'kt = ' + str(kt)
    raw_input("press enter to continue")
    #mlab.show()
    #if (kt == 0):
    #    raw_input("press enter to continue")
        
