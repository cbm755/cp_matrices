from itertools import izip
from math import sqrt
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree
#from parallel_kdquery import cKDTree_MP
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab

from cp.tools.io import load_ply
from cp.surfaces.triangulation_fast import FindClosestPointToTriSet


class Mesh(object):
    """Mesh that knows how to compute its closest point representation."""
    def __init__(self, vertices, faces, nprocs=8, parallel=False):
        self.vertices = vertices
        self.faces = faces
        self.nprocs = nprocs
        if parallel:
            self.tree = cKDTree_MP(self.vertices, self.nprocs)
        else:
            self.tree = cKDTree(self.vertices)
        self.vertices2faces = build_vertices2faces(self.faces)

    def grid(self, num_blocks=100, levels=2):
        #bounding_box = [self.vertices.min(axis=0), self.vertices.max(axis=0)]
        grid, initial_dx = build_grid(num_blocks,
                                      3*[float(int(2*self.vertices.min()))],
                                      3*[float(int(2*self.vertices.max()))])
        grid = grid.reshape((3, -1)).T
        final_dx = initial_dx / 3.**(levels-1)
        final_bandwidth = np.max(bandwidth(final_dx, dim=3))
        level_dx = initial_dx.copy()
        if levels > 1:
            level_bandwidth = sqrt(np.dot(level_dx,
                                          level_dx)) / 2. + np.max(bandwidth(level_dx, dim=3)) # final_bandwidth
        else:
            level_bandwidth = final_bandwidth
        offsets = np.mgrid[-1:1:3j,
                           -1:1:3j,
                           -1:1:3j].reshape((3, -1))
        for i in xrange(levels):
            print "level_dx", level_dx
            print "level_bandwidth", level_bandwidth
            index, grid, dist = self._new_level(grid, level_dx, level_bandwidth)

            # Last loop
            if i >= (levels-1) or levels == 1:
                pass  # Exits the loop, we're done
            else:
                level_dx = initial_dx / 3.**(i+1)
                # Last but one loop (the case levels==1 is handled by
                # the if above the loop)
                if i == (levels-2):  
                    level_bandwidth = final_bandwidth
                # Next loop is not the final one: different bandwidth
                else:
                    level_bandwidth = (sqrt(np.dot(level_dx, level_dx)) / 2. +
                                       np.max(bandwidth(level_dx, dim=3)))
                # Refine the grid
                grid = (grid[..., np.newaxis] +
                        (offsets * level_dx[:, np.newaxis])[np.newaxis, ...])
                grid = grid.transpose((0, 2, 1)).reshape((-1, 3))
                
        return index, grid, dist, final_dx

    def _new_level(self, grid, level_dx, level_bandwidth):
        # Giving it the upper bound makes it *much* faster
        dist, index = self.tree.query(grid, distance_upper_bound=level_bandwidth)
        is_within_bandwidth = np.isfinite(dist)
        index = index[is_within_bandwidth]
        grid = grid[is_within_bandwidth]
        dist = dist[is_within_bandwidth]
        return index, grid, dist

    def _refine(self, index, grid):
        return refine(index, grid, self.vertices2faces,
                      self.vertices, self.faces)
    
    def closest_point(self, index, grid):
        res = self._refine(index, grid)
        cp = res[:, 1:]
        dist = np.sqrt(res[:, 0])
        return cp, dist, 0, {}

    def plot(self):
        mlab.triangular_mesh(self.vertices[:,0],
                             self.vertices[:,1],
                             self.vertices[:,2],
                             self.faces, opacity=0.3)

    def plot_grid(self, grid, half=True, colour=(0.2, 0.8, 0.8)):
        if half:
            mask = grid[:,0] > 0
        else:
            mask = slice(None)
        mlab.points3d(grid[mask,0], grid[mask,1], grid[mask,2], 
                      mode='point', color=colour)


def build_grid(n, ll, ur, filled=True):
    """
    Makes a multi-dimensional coarse grid.

    It does not special-case the 2D case.

    Input
    -----
    n : scalar
        number of blocks in each dimension
    ll : array-like
         "lowest, left most" corner of the mesh
    ur : array-like
         "upper most, right most" corner of the mesh
    filled : bool, default is True
             whether to return a filled or an open grid

    Output
    ------
    grid : array
           desired grid
    dx : array
         grid spacing
    """
    ll, ur = np.asarray(ll, dtype=np.float64), np.asarray(ur, dtype=np.float64)
    f = np.mgrid if filled else np.ogrid
    grid = f[tuple(slice(ll_i, ur_i, complex(0, n))
                   for ll_i, ur_i in zip(ll, ur))]
    dx = (ur - ll) / (n-1)
    return grid, dx

def bandwidth(dx, dim, p=3, diff_stencil_arm=1):
    """Bandwith value.

    Input
    -----
    dx : scalar or array
         grid spacing
    p : scalar
        interpolation order
    diff_stencil_arm : scalar, default is 1 (second order centered
                       difference Laplacian)
                       lenght of the differencing stencil arm
        
    Reference
    ---------
    Ruuth & Merriman
    """
    # 1.0001 is a security factor
    lambda_ = 1.0001 * sqrt((dim-1.0)*((p+1)/2.0)**2 + 
                        (diff_stencil_arm + (p+1)/2.0)**2) * dx
    return lambda_

def build_vertices2faces(faces):
    """Dict that maps each vertex to the indices of all its faces."""
    vertices2faces = defaultdict(list)
    for i, face in enumerate(faces):
        for vertex in face:
            vertices2faces[vertex].append(i)
    return vertices2faces

def refine(index_cp, grid_fine, vertex2faces, vertices, faces):
    """Actual closest points to a triangulated surface.

    Given the closest vertex to any point, finds the actual closest
    point in the triangulated surface.

    For each given point, whose closest point we want to find, the
    algorithm finds the closest point to all the faces adjacent to the
    closest vertex, and then selects the minimum.

    Input
    -----
    index_cp : array
               Contains the indices of the closest vertices
    grid_fine : array
                Points whose closest points we want to find
    vertex2faces : dict
                   Given by build_vertices2faces
    vertices, faces : arrays, given by load_ply

    Output
    ------
    cp : array of shape (npoints, 4)
         Coordinates of the closest points, and squareddistance from
         its original point.
         cp[:, 0] is the squared distance (redundant information)
         cp[:, 1:] are the coordinates of the closest points
    """
    N = len(index_cp)
    cp = np.empty((N, 4), dtype=np.float64)
    for i, (index_of_cp, point) in enumerate(izip(index_cp, grid_fine)):
        face_indices_that_share_that_vertex = vertex2faces[index_of_cp]
        # Much faster using np.take than fancy indexing
        possible_faces = np.take(faces,
                                 face_indices_that_share_that_vertex,
                                 axis=0)
        # Cythonize this and it gets way faster
        cp[i] = FindClosestPointToTriSet(point[0], point[1], point[2],
                                         possible_faces,
                                         vertices)
    return cp

if __name__ == '__main__':
    v, f = load_ply('cp/tests/data/Armadillo.ply')
    m = Mesh(v, f)

    v2, f2 = load_ply('cp/tests/data/eight_refined.ply')
    m2 = Mesh(v2, f2)
    index, grid, distance, dx = m2.grid(num_blocks=10, levels=4)
    cp, dist, _, _ = m2.closest_point(index, grid)
