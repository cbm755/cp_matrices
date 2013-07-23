function [pass, str] = test_refine_cpgrid_bw_nd()
  str = ['refinement test: ndgrid bw refinements: regression tests'];
  pass = [];
  c = 0;

  dx = 1;

  % make vectors of x, y, positions of the grid
  x = (-3.0:dx:3.0)';
  dim = 4;
  x1d = {};
  for d = 1:dim
    x1d{d} = x;
  end

  xx = {};
  [xx{1:dim}] = ndgrid(x1d{1:dim});
  for d = 1:dim
    xx{d} = xx{d}(:);
  end

  cen = dx/3*ones(1,dim);
  cpfun = @(x) cpCircleInHighDim(x, 1, cen);

  g0.dim = dim;
  g0.dx = dx;
  g0.x1d = x1d;
  g0.cpfun = cpfun;
  g0.band = 1:length(xx{1});
  g0.x = xx;
  [g0.cpx,g0.dist] = cpfun(xx);

  bw = rm_bandwidth(dim, 3, 1, 1.0001);

  % need to refine a bit with stencil code until it fits in the virtual
  % ndgrid, need a wide stencil for the effective bw g1 to be big
  % enough
  g1 = g0;
  g1 = refine_cpgrid(g1, 3);
  g1 = refine_cpgrid(g1, 6);
  g2 = refine_cpgrid_bw(g1, bw);
  g3 = refine_cpgrid_bw(g2, bw);
  g4 = refine_cpgrid_bw(g3, bw);

  c = c + 1;
  pass(c) = length(g2.band) == 20186;
  c = c + 1;
  pass(c) = length(g3.band) == 40430;
  c = c + 1;
  pass(c) = length(g4.band) == 80816;

  g5 = g0;
  g5 = refine_cpgrid(g5, 1);
  g5 = refine_cpgrid(g5, 1);
  g5 = refine_cpgrid(g5, 5);
  g5 = refine_cpgrid_bw(g5, bw);
  g5 = refine_cpgrid_bw(g5, bw);

  c = c + 1;
  pass(c) = length(g5.band) == length(g4.band);

