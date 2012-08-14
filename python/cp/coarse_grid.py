from __future__ import division

from math import sqrt
import numpy as np



class CoarseGrid(object):
    """Coarse grid of a surface that knows which blocks are close
    enough to the surface that they need a finer mesh"""
    def __init__(self, surface):
        self.surface = surface
        
    def grid(self, n, ll, ur, filled=True):
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
        """
        ll, ur = np.asarray(ll, dtype=np.float64), np.asarray(ur, dtype=np.float64)
        f = np.mgrid if filled else np.ogrid
        grid = f[tuple(slice(ll_i, ur_i, complex(0, n))
                       for ll_i, ur_i in zip(ll, ur))]
        self.dx = (ur - ll) / (n-1)
        return grid

    def bandwidth(self, p, diff_stencil_arm=1):
        """Bandwith value.

        Input
        -----
        p : scalar
            interpolation order
        diff_stencil_arm : scalar, default is 1 (second order centered
                           difference Laplacian)
                           lenght of the differencing stencil arm
            
        Reference
        ---------
        Ruuth & Merriman
        """
        dim = self.surface._dim
        dx = self.dx  # TODO fine dx
        # 1.0001 is a security factor
        lam = 1.0001 * sqrt((dim-1.0)*((p+1)/2.0)**2 + 
                            (diff_stencil_arm + (p+1)/2.0)**2) * dx
        return lam

    def bound(self, bandwidth):
        """'Half a diagonal' of a hypercube + bandwidth

        If the distance to the closest point from the center point in
        a hypercube is further than this bound, we don't need to
        refine that hypercube."""
        dx = self.dx
        return sqrt(np.dot(dx, dx)) / 2. + np.max(bandwidth)

    def mask(self, grid, bandwidth):
        _, dist, _, _ = self.surface.closestPointToCartesianGrid(grid)
        return np.abs(dist) <= self.bound(bandwidth)

    def build_index_mappers(self, mask):
        self.linear_to_grid = {i:j for i, j in enumerate(zip(*np.nonzero(mask)))}
        self.grid_to_linear = {v:k for k, v in self.linear_to_grid.iteritems()}

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from mayavi import mlab
    
    from surfaces import Circle, Hemisphere

    circle = True
    if circle:
        sur = Circle()
        c = CoarseGrid(sur)
        g = c.grid(30., [-2., -2.], [2., 2.])
        b = c.bandwidth(3., 1.)
        mask = c.mask(g, b)
        c.build_index_mappers(mask)
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        scatter = ax.scatter(g[0], g[1], c=mask, cmap=cm.flag, s=10, linewidth=0)
        heat = ax.imshow(mask, extent=[-2, 2, -2, 2], interpolation='nearest')
        x, y = sur.ParamGrid()
        ax.scatter(x, y, s=2)
        plt.show()

    else:
        sur = Hemisphere()
        c = CoarseGrid(sur)
        g = c.grid(40, [-2, -2, -2], [2, 2, 2])
        b = c.bandwidth(3., 1.)
        mask = c.mask(g, b)
        c.build_index_mappers(mask)
        fig = mlab.figure()
        scalars = np.ones_like(mask, dtype=np.int)
        scalars[mask] = 2.
        s = mlab.points3d(g[0].ravel(), g[1].ravel(), g[2].ravel(), scalars.ravel(), scale_factor=0.01, transparent=True)
        s = mlab.contour3d(g[0], g[1], g[2], scalars, opacity=0.7, contours=2, color=(0.2, 0.2, 0.2))
        x, y, z = sur.ParamGrid()
        s2 = mlab.mesh(x, y, z, color=(0.2, 0.8, 0.2))
        mlab.show()
