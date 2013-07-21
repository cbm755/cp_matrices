function gg = refine_cpgrid_sten_2d(g, p)
%REFINE_CPGRID_STEN_2D CP grid refinement in 2D using interp stencil
%
%   Private: you probably want refine_grid()
%
%   gg = refine_cpgrid_sten_2d(g, p)
%     uses (stencils for) degree p interpolation.  p=3 is the
%     default if omitted.
%
%   TODO: doesn't deal with open surfaces properly (bdy)
%   TODO: could be optimized quite a bit if computing
%         CPs is expensive
%   TODO: 

  if nargin < 2
    p = 3;
  end
  assert(g.dim == 2);

  cpfun = g.cpfun;
  dim = g.dim;

  % build a new grid with half the grid spacing of the old one
  dx = g.dx / 2;
  relpt = [g.x1d(1) g.y1d(1)];
  x1d  = ( g.x1d(1):dx:g.x1d(end) )';
  y1d  = ( g.y1d(1):dx:g.y1d(end) )';
  nx = length(x1d);
  ny = length(y1d);


  %% The main idea
  % build an interpolation matrix in the new grid for the old
  % closest points.  If dx hasn't changed much (and p large
  % enough), the column space will be the band.
  %
  % todo: is it worth having a separate routine?  Not much faster,
  % maybe 10% or 20%, see nd version
  [tilde,Ej,tilde] = interp2_matrix(x1d, y1d, g.cpx, g.cpy, p);
  band = unique(Ej);

  % TODO: cpgrid objects should carry these conversion functions around
  [J,I] = ind2sub([ny nx], band);

  x = relpt(1) + (I-1)*dx;
  y = relpt(2) + (J-1)*dx;

  % find the new closest points.  Some waste here b/c we probably
  % knew some of them already.
  [cpx, cpy, dist] = cpfun(x, y);


  %% Safety check
  % if we're working with minimal bands with no safety margin, its
  % possible the above could miss a few grid points (I think).
  % We can do it again using the new closest points.
  % TODO: perhaps needs a toggle
  [tilde,Ej,tilde] = interp2_matrix(x1d, y1d, cpx, cpy, p);
  band2 = unique(Ej);

  if length(band) == length(band2) && all(band == band2)
    % gave same thing, no need to recompute cp's.
  else
    % TODO: could merge instead of recomputing all those CPs
    band = band2;
    [J,I] = ind2sub([ny nx], band);
    x = relpt(1) + (I-1)*dx;
    y = relpt(2) + (J-1)*dx;
    [cpx, cpy, dist] = cpfun(x, y);
  end

  gg.dim = dim;
  gg.dx = dx;
  gg.x1d = x1d;
  gg.y1d = y1d;
  gg.cpfun = cpfun;
  gg.band = band;
  gg.x = x;
  gg.y = y;
  gg.cpx = cpx;
  gg.cpy = cpy;
  gg.dist = dist;
