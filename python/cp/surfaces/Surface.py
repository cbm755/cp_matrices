from math import sqrt

import numpy as np


class Surface(object):
    def closest_point(self, x):
        raise NotImplementedError("Should be implemented in subclass")
    
    def grid(self, num_blocks_per_dim=100, levels=2, p=3, diff_stencil_arm=1):
        grid, initial_dx = build_grid(num_blocks_per_dim,
                                      2 * self.bounding_box[0],
                                      2 * self.bounding_box[1])
        grid = grid.reshape((self.dim, -1)).T
        final_dx = initial_dx / 3.**(levels-1)
        final_bandwidth = np.max(bandwidth(final_dx,
                                           p=p,
                                           diff_stencil_arm=diff_stencil_arm,
                                           dim=self.dim))
        level_dx = initial_dx.copy()
        if levels > 1:
            level_bandwidth = sqrt(np.dot(level_dx, level_dx)) / 2. +\
              np.max(bandwidth(level_dx,
                               p=p,
                               diff_stencil_arm=diff_stencil_arm,
                               dim=self.dim))
        else:
            level_bandwidth = final_bandwidth
        offsets = np.mgrid[tuple(
            slice(-1, 1, complex(0, 3)) for _ in xrange(self.dim)
            )].reshape((self.dim, -1))
        for i in xrange(levels):
            print "level_dx", level_dx
            print "level_bandwidth", level_bandwidth
            # index if triangulated surface, else cp
            index, dist, grid = self._new_level(grid, level_dx, level_bandwidth)
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
                                       np.max(bandwidth(level_dx,
                                                        p=p,
                                                        diff_stencil_arm=diff_stencil_arm,
                                                        dim=self.dim)))
                # Refine the grid
                grid = (grid[..., np.newaxis] +
                        (offsets * level_dx[:, np.newaxis])[np.newaxis, ...])
                grid = grid.transpose((0, 2, 1)).reshape((-1, self.dim))

        # index if triangulated surface, else cp
        return index, dist, grid, final_dx

    def _new_level(self, grid, level_dx, level_bandwidth):
        cp, dist, _, _ = self.closest_point(grid)
        is_within_bandwidth = dist <= level_bandwidth
        cp = cp[is_within_bandwidth]
        grid = grid[is_within_bandwidth]
        dist = dist[is_within_bandwidth]
        return cp, dist, grid

class ShapeWithBdy(Surface):
    pass

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
