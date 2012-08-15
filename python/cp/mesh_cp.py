from collections import defaultdict
from itertools import izip

import numpy as np
import scipy.spatial as ss

from cp.surfaces.triangulation_fast import FindClosestPointToTriSet as fast_FindClosestPointToTriSet


def build_vertices2faces(faces):
    """Returns dict that maps each vertex to the indexes of all its faces."""
    vertices2faces = defaultdict(list)
    for i, face in enumerate(faces):
        for vertex in face:
            vertices2faces[vertex].append(i)
    return vertices2faces

def refine(index_cp, grid_fine, vertex2faces, vertices, faces):
    N = len(index_cp)
    cp = np.empty((N, 4), dtype=np.float64)
    for i, (index_of_cp, point) in enumerate(izip(index_cp, grid_fine)):
        faces_indexes_that_share_that_vertex = vertex2faces[index_of_cp]
        #possible_faces = faces[faces_indexes_that_share_that_vertex]
        # Much faster using np.take
        possible_faces = np.take(faces,
                                 faces_indexes_that_share_that_vertex,
                                 axis=0)
        # Cythonize this and it gets 15x faster
        cp[i] = fast_FindClosestPointToTriSet(point[0], point[1], point[2],
                                              possible_faces,
                                              vertices)
    return cp


class MeshCP(object):
    def __init__(self, vertices, faces):
        self.vertices = vertices
        self.faces = faces
        self.tree = ss.cKDTree(self.vertices)
        self.vertices2faces = build_vertices2faces(self.faces)
        
    def query(self, points, bound=np.inf):
        dist, index = self.tree.query(points, distance_upper_bound=bound)
        cp = refine(index, points, self.vertices2faces,
                    self.vertices, self.faces)
        return cp  # Returns (square distance, cp_x, cp_y, cp_z)
