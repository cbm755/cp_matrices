from cp.coarse_grid import CoarseGrid
from cp.surfaces import Circle

class test_CoarseGrid():
    def test_index_mappings(self):
        c = CoarseGrid(Circle())
        g = c.grid(30., [-2., -2.], [2., 2.])
        b = c.bandwidth(3., 1.)
        mask = c.mask(g, b)
        c.build_index_mappers(mask)
        for k in c.linear_to_grid:
            assert k == c.grid_to_linear[c.linear_to_grid[k]]
