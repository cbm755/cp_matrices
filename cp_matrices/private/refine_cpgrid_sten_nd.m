function gg = refine_cpgrid_sten_nd(g, p)
%REFINE_CPGRID_STEN_ND CP grid refinement in n-dimensions using stencil
%
%   Private: you probably want refine_grid()

  % see refine_cpgrid_sten_2d for comments

  if nargin < 2
    p = 3;
  end

  cpfun = cp.cpfun;
  dim = cp.dim;
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
  %[Ei,Ej,Es] = interpn_matrix(cp2.x1d, cp.cpx, p);
  %band2 = unique(Ej);
  %toc
  tic
  band = findn_band(cp2.x1d, cp.cpx, p);
  toc
  %assert(all(band == band2))

  gg.band = band;

  [I{1:dim}] = ind2sub(Ns, band);

  % build the xyz coordinates of the points in the band
  for d = 1:dim
    x{d} = relpt(d) + (I{d}-1)*dx;
  end

  % find the closest point
  [cpx,dist] = cpfun(gg.x);

  gg.x = x;
  gg.cpx = cpx;
  gg.dist = gdist;

