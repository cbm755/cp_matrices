from itertools import izip
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree
#from parallel_kdquery import cKDTree_MP
#try:
#    from mayavi import mlab
#except ImportError:
#    from enthought.mayavi import mlab

from cp.tools.io import load_ply
from cp.surfaces.triangulation_fast import FindClosestPointToTriSet

from Surface import Surface


class Mesh(Surface):
    """Mesh that knows how to compute its closest point representation."""
    def __init__(self, vertices, faces, nprocs=8, parallel=False):
        self.vertices = vertices
        self.faces = faces
        self.nprocs = nprocs
        self.dim = 3
        self.bounding_box = [np.array(self.dim * [vertices.min()]),
                             np.array(self.dim * [vertices.max()])]
        if parallel:
            self.tree = cKDTree_MP(self.vertices, self.nprocs)
        else:
            self.tree = cKDTree(self.vertices)
        self.vertices2faces = build_vertices2faces(self.faces)

    def _new_level(self, grid, level_dx, level_bandwidth):
        # Giving it the upper bound makes it *much* faster
        dist, index = self.tree.query(grid, distance_upper_bound=level_bandwidth)
        is_within_bandwidth = np.isfinite(dist)
        index = index[is_within_bandwidth]
        grid = grid[is_within_bandwidth]
        dist = dist[is_within_bandwidth]
        return index, dist, grid

    def _refine(self, grid, index):
        res = refine(grid, index, self.vertices2faces,
                     self.vertices, self.faces)
        cp = res[:, 1:]
        dist = np.sqrt(res[:, 0])
        return cp, dist
    
    def closest_point(self, grid, index):
        """index is returned by self.grid. It contains information
        about the closest vertices.
        """
        cp, dist = self._refine(grid, index)
        return cp, dist, 0, {}

def build_vertices2faces(faces):
    """Dict that maps each vertex to the indices of all its faces."""
    vertices2faces = defaultdict(list)
    for i, face in enumerate(faces):
        for vertex in face:
            vertices2faces[vertex].append(i)
    return vertices2faces

def refine(grid_fine, index_cp, vertex2faces, vertices, faces):
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
    index, distance, grid, dx = m2.grid(num_blocks=10, levels=4)
    cp, dist, _, _ = m2.closest_point(index, grid)
