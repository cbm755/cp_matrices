function gg = refine_cpgrid_sten_2d(g, p)
%REFINE_CPGRID_STEN_2D CP grid refinement in 2D using interp stencil
%
%   Private: you probably want refine_cpgrid()
%
%   gg = refine_cpgrid_sten_2d(g, p)
%     uses (stencils for) degree p interpolation.  p=3 is the
%     default if omitted.
%
%   TODO: doesn't deal with open surfaces properly (bdy)

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


  %% Iteration
  % if we're working with minimal bands with no safety margin, its
  % possible the above could miss a few grid points.  We repeat the
  % above as an iteration and look for a fixed point (usually one
  % step).
  iter = 1;
  while(1)
    [tilde,Ej,tilde] = interp2_matrix(x1d, y1d, cpx, cpy, p);
    band2 = unique(Ej);

    if length(band) == length(band2) && all(band == band2)
      fprintf('  refine iter%d: no change, stopping iteration\n', iter);
      break
    else
      %% Compute new closest points
      % this is simple but needs to compute many CPs, instead we only
      % compute the new ones
      %band = band2;
      %[J,I] = ind2sub([ny nx], band);
      %x = relpt(1) + (I-1)*dx;
      %y = relpt(2) + (J-1)*dx;
      %[cpx, cpy, dist] = cpfun(x, y);

      [tilde,ii] = setdiff(band2,band);
      [tilde,i4,iii] = intersect(band,band2);
      [tilde,i5] = setdiff(band,band2);
      if isempty(i5)
        fprintf('  refine iter%d: found %d new points\n', iter, length(ii));
      else
        fprintf('  refine iter%d: found %d new points (and lost %d others)\n', ...
                iter, length(ii), length(i5));
      end
      [J,I] = ind2sub([ny nx], band2);
      [J2,I2] = ind2sub([ny nx], band2(ii));
      % build the xyz coordinates of the points in the band
      x = relpt(1) + (I-1)*dx;
      y = relpt(2) + (J-1)*dx;
      xT = relpt(1) + (I2-1)*dx;
      yT = relpt(2) + (J2-1)*dx;
      % find the closest point
      [cpxT, cpyT, distT] = cpfun(xT, yT);
      dist2 = zeros(size(band2));
      dist2(iii) = dist(i4);
      dist2(ii) = distT;
      cpx2 = zeros(size(band2));
      cpx2(iii) = cpx(i4);
      cpx2(ii) = cpxT;
      cpy2 = zeros(size(band2));
      cpy2(iii) = cpy(i4);
      cpy2(ii) = cpyT;
      if (1==0)  % debugging, check it matches recomputing everything
        [cpx3, cpy3, dist3] = cpfun(x, y);
        assert(norm(dist3-dist2)==0)
        assert(norm(cpx3-cpx2)==0)
        assert(norm(cpy3-cpy2)==0)
      end
      band = band2;
      cpx = cpx2;
      cpy = cpy2;
      dist = dist2;
    end
    iter = iter + 1;
  end % while

  fprintf('  refine: found %d points in %d iterations\n', length(band), iter);
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
