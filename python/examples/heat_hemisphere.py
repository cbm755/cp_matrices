r"""
Solve the diffusion (heat) equation

     $u_t = \kappa \triangle_s u$

on a surface.  $\triangle_s$ is the Laplace--Beltrami operator.
$\kappa$ is a scalar coefficient.

Uses implicit backward Euler timestepping.
"""

import numpy as np
from scipy.sparse import identity
from scipy.sparse.linalg import gmres, spsolve
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab
from time import time

from cp import surfaces, cpGrid, cpOps


class HeatEquation(object):
    def __init__(self, surface, has_bdy):
        self.surface = surfaces.CPBar(surface) if has_bdy else surface
        self.has_bdy = has_bdy

    def grid(self, x0, initial_dx, max_level):
        # base grid point
        self.x0 = x0
        # Grid spacing from which it'll start refining
        self.initial_dx = initial_dx
        # Maximum levels of refinement
        self.max_level = max_level
        
        TreeGrid = cpGrid.CPGrid('Spamname', self.surface.cp,
                                 self.surface._dim, self.x0,
                                 self.initial_dx,
                                 self.interp_degree,
                                 levels=self.max_level+1)

        # Hemisphere parametrization in cartesian coordinates, for example for
        # plotting
        self.x, self.y, self.z = surface.ParamGrid(64)
        self.PlotPts = np.vstack((x.ravel(), y.ravel(), z.ravel())).T

        TreeGrid.findStencilSets(self.max_level)
        TreeGrid.buildListsFromStencilSets(self.max_level)

        def boundary_function(bdy):
            if bdy == 1 or bdy == 2:
                return 'dirichlet_2nd_order'
            else:
                return None

        TreeGrid.findStencilsOnStencilSets(self.max_level,
                                           boundary_function)
        

        self.Lev = TreeGrid.Levolve[self.max_level]
        self.dx = self.Lev[0].dx
        self.Lex = TreeGrid.Lextend[self.max_level]
        self.Grid = TreeGrid.Grids[self.max_level]

    def buildDiffMatrix(self):
        return cpOps.buildDiffMatrix(self.Lev, self.Lex)
    
    def buildExtensionMatrix(self):
        return cpOps.buildExtensionMatrix(self.Lev, self.Lex)

    def buildEPlotMatrix(self):
        return cpOps.buildEPlotMatrix(self.Grid, self.Lev, self.Lex, 
                                            self.PlotPts, self.interp_degree)

    def buildLinearDiagonalSplitting(self, D, E):
        return cpOps.LinearDiagonalSplitting(D, E)

    def solve(self, implicit, Nsteps, dt, kappa, u0):
        D = buildDiffMatrix()
        E = buildExtensionMatrix()
        EPlot = buildEPlotMatrix()
        if implicit:
            M = buildLinearDiagonalSplitting(D, E)
        
    
if __name__ == '__main__':
    h = surfaces.Hemisphere()
    ex = HeatEquation(h, has_bdy=True)
    ex.grid(x0=np.zeros(h._dim), initial_dx=4, max_level=4)
    ex.solve(True, 100, 0.1 * ex.dx, 1, np.random.randn(
