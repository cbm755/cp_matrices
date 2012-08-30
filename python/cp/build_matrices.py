from operator import mul
import numpy as np
from scipy.sparse import coo_matrix

from cp.cpOps import findGridInterpBasePt, buildInterpWeights


def build_diff_matrix(int_band, dx, shape):
    dim = len(shape)
    if dim == 3:
        # Assume now that band is not a linear index, but a 3D index (ie, int_band)
        weights = np.array([-6, 1, 1, 1, 1, 1, 1]) / np.max(dx)**2
        offsets = np.array([[0,0,0],
                            [0,0,1],
                            [0,1,0],
                            [1,0,0],
                            [0,0,-1],
                            [0,-1,0],
                            [-1,0,0]])
        Li = np.tile(np.arange(int_band.shape[0])[:, np.newaxis], weights.size)
        Lj = np.zeros_like(Li)
        Ls = np.zeros_like(Li, dtype=np.float)
        i, j, k = int_band.T  #np.ravel_multi_index(band.T, shape)

        for c in xrange(weights.size):
            ii = i + offsets[c, 0]
            jj = j + offsets[c, 1]
            kk = k + offsets[c, 2]
            Ls[:, c] = weights[c]
            Lj[:, c] = np.ravel_multi_index((ii, jj, kk), shape)
        L = coo_matrix((Ls.ravel(), (Li.ravel(), Lj.ravel())),
                       (int_band.shape[0], shape[0] * shape[1] * shape[2])).tocsr()
        return L[:, np.ravel_multi_index(int_band.T, shape)]
    elif dim == 2:
        weights = np.array([-4, 1, 1, 1, 1]) / np.max(dx)**2
        offsets = np.array([[0, 0],
                            [0, 1],
                            [1, 0],
                            [0, -1],
                            [-1, 0]])
        Li = np.tile(np.arange(int_band.shape[0])[:, np.newaxis], weights.size)
        Lj = np.zeros_like(Li)
        Ls = np.zeros_like(Li, dtype=np.float)
        i, j = int_band.T  #np.ravel_multi_index(band.T, shape)

        for c in xrange(weights.size):
            ii = i + offsets[c, 0]
            jj = j + offsets[c, 1]
            Ls[:, c] = weights[c]
            Lj[:, c] = np.ravel_multi_index((ii, jj), shape)
        L = coo_matrix((Ls.ravel(), (Li.ravel(), Lj.ravel())),
                       (int_band.shape[0], shape[0] * shape[1])).tocsr()
        return L[:, np.ravel_multi_index(int_band.T, shape)]

def build_interp_matrix(int_band, cp, dx, p, ll, shape):
    dim = len(shape)
    offsets = np.mgrid[tuple(
        slice(p+1) for _ in xrange(dim))].reshape((dim, -1))
    Ei = np.tile(np.arange(cp.shape[0])[:, np.newaxis], (p+1)**dim)
    Ej = np.zeros_like(Ei)
    Es = np.zeros_like(Ei, dtype=np.float)

    base_points = findGridInterpBasePt(cp, dx, ll, p)
    weights = buildInterpWeights(base_points * dx + ll, cp, dx, p+1)
    Ej = np.ravel_multi_index((base_points[..., np.newaxis] +
                               offsets[np.newaxis, ...]).transpose(
                                   (0, 2, 1)).reshape((-1, dim)).T,
                                   shape)
    Es = weights
    E = coo_matrix((Es.ravel(), (Ei.ravel(), Ej.ravel())),
                   (cp.shape[0], reduce(mul, shape))).tocsr()
    return E[:, np.ravel_multi_index(int_band.T, shape)]
