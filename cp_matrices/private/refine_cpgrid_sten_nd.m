function gg = refine_cpgrid_sten_nd(g, p)
%REFINE_CPGRID_STEN_ND CP grid refinement in n-dimensions using stencil
%
%   Private: you probably want refine_cpgrid()

  % see refine_cpgrid_sten_2d for comments

  if nargin < 2
    p = 3;
  end

  cpfun = g.cpfun;
  dim = g.dim;
  assert(dim == length(g.x1d));

  relpt = zeros(1, dim);
  for d = 1:dim
    relpt(d) = g.x1d{d}(1);
  end

  dx = g.dx / 2;

  Ns = zeros(1, dim);
  for d = 1:dim
    x1d{d} = (g.x1d{d}(1):dx:g.x1d{d}(end))';
    Ns(d) = length(x1d{d});
  end
  gg.dim = dim;
  gg.dx = dx;
  gg.x1d = x1d;
  gg.cpfun = cpfun;


  %tic
  %[Ei,Ej,Es] = interpn_matrix(gg.x1d, g.cpx, p);
  %band2 = unique(Ej);
  %toc
  tic
  band = findn_band(gg.x1d, g.cpx, p);
  toc
  %assert(all(band == band2))

  [I{1:dim}] = ind2sub(Ns, band);

  % build the xyz coordinates of the points in the band
  x = {};
  for d = 1:dim
    x{d} = relpt(d) + (I{d}-1)*dx;
  end

  % find the closest point
  [cpx,dist] = cpfun(x);

  iter = 1;
  while(1)
    band2 = findn_band(gg.x1d, cpx, p);

    if length(band) == length(band2) && all(band == band2)
      fprintf('  refine iter%d: no change, stopping iteration\n', iter);
      break
    else
      [tilde,ii] = setdiff(band2,band);
      [tilde,i4,iii] = intersect(band,band2);
      [tilde,i5] = setdiff(band,band2);
      if isempty(i5)
        fprintf('  refine iter%d: found %d new points\n', iter, length(ii));
      else
        fprintf('  refine iter%d: found %d new points (and lost %d others)\n', ...
                iter, length(ii), length(i5));
      end
      [I{1:dim}] = ind2sub(Ns, band2);
      [I2{1:dim}] = ind2sub(Ns, band2(ii));
      % build the xyz coordinates of the points in the band
      x = {};
      xT = {};
      for d = 1:dim
        x{d} = relpt(d) + (I{d}-1)*dx;
        xT{d} = relpt(d) + (I2{d}-1)*dx;
      end
      % find the closest point
      [cpxT, distT] = cpfun(xT);
      dist2 = zeros(size(band2));
      dist2(iii) = dist(i4);
      dist2(ii) = distT;
      cpx2 = {};
      for d = 1:dim
        cpx2{d} = zeros(size(band2));
        cpx2{d}(iii) = cpx{d}(i4);
        cpx2{d}(ii) = cpxT{d};
      end

      if (1==0)  % debugging, check it matches recomputing everything
        [cpx3, dist3] = cpfun(x);
        assert(norm(dist3-dist2)==0)
        for d = 1:dim
          assert(norm(cpx3{d}-cpx2{d})==0)
        end
      end
      band = band2;
      cpx = cpx2;
      dist = dist2;
    end
    iter = iter + 1;
  end % while

  fprintf('  refine: found %d points in %d iterations\n', length(band), iter);
  gg.band = band;
  gg.x = x;
  gg.cpx = cpx;
  gg.dist = dist;

