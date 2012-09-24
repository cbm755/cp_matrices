
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


PLOT = True

# Load vertices and faces, and instantiate surface
vert, faces = load_ply('brain-lh_scale_1.ply')
#m = Mesh(vert, faces)

if PLOT:
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


(dx, Tf, numtimesteps, dt, initial_u_plot) = pickle.load(file('brain_IC_plotvec.pickle'))
(dx, Tf, numtimesteps, dt, uplot) = pickle.load(file('brain_final_plotvec.pickle'))


if PLOT:
    src.data.point_data.scalars = uplot
    src.data.point_data.scalars.name = 'scalars'
    src.data.modified()
