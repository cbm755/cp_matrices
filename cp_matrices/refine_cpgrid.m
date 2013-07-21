function gg = refine_cpgrid(g, p)
%REFINE_CPGRID  CP grid refinement
%
%   cpgrid_new = refine_cpgrid(cpgrid_old)
%
%   cpgrid_new = refine_cpgrid(cpgrid_old, p)
%     uses (stencils for) degree p interpolation.  p=3 is the
%     default if omitted.
%
%   cpgrid_new = refine_cpgrid(cpgrid_old, 'bw', bw)
%     uses a bandwidth of bw*dx instead.
%     TODO: not implemented, although see "refine_grid", an earlier
%     code that does this.

  if nargin < 2
    p = 3;
  end

  if iscell(g.x1d)
    gg = refine_cpgrid_nd_degp(g, p);
  elseif g.dim == 2
    gg = refine_cpgrid_2d_degp(g, p);
  elseif g.dim == 3
    gg = refine_cpgrid_3d_degp(g, p);
  else
    error('this grid object cannot be handled');
  end

