function gg = refine_cpgrid(g, p)
%REFINE_CPGRID  CP grid refinement using stencils
%
%   cpgrid_new = refine_cpgrid(cpgrid_old)
%   cpgrid_new = refine_cpgrid(cpgrid_old, p)
%     Given an existing, possible banded grid, construct a new one
%     with dx/2.  The new band is large enough for the stencils of
%     degree p interpolation (p=3 by default if omitted).  The band
%     consists of points in the column-space of the corresponding
%     degree p interpolation operator.
%
%   cpgrid_new = refine_cpgrid(cpgrid_old, 'bw', bw)
%     uses a bandwidth of bw*dx instead.
%     TODO: not implemented, although see "refine_cpgrid_bw", for a
%     refinment code which does this.
%
%   Tips: if you want to refine more than once, it should be safe
%   to use a lower p for the intermediate grids.  Something like this:
%     g = refine_cpgrid(g, 1)
%     g = refine_cpgrid(g, 1)
%     g = refine_cpgrid(g, 1)
%     g = refine_cpgrid(g, 5)

  if nargin < 2
    p = 3;
  end

  if iscell(g.x1d)
    gg = refine_cpgrid_sten_nd(g, p);
  elseif g.dim == 2
    gg = refine_cpgrid_sten_2d(g, p);
  elseif g.dim == 3
    gg = refine_cpgrid_sten_3d(g, p);
  else
    error('this grid object cannot be handled');
  end

