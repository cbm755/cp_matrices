function [pass, str] = test_refine_cpgrid_nd()
  str = ['refinement test: ndgrid stencil refinements: regression tests'];
  pass = [];
  c = 0;
  
  dx = 1;

  % make vectors of x, y, positions of the grid
  x = (-2.0:dx:2.0)';
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
  g0.cpx = cpfun(xx);
  g0.x = xx;

  c = c + 1;
  pass(c) = length(g0.band) == 625;


  g1 = g0;
  g1 = refine_cpgrid(g1, 3);

  c = c + 1;
  pass(c) = length(g1.band) == 976;

  g1 = refine_cpgrid(g1, 3);
  g1 = refine_cpgrid(g1, 3);
  g1 = refine_cpgrid(g1, 3);

  c = c + 1;
  pass(c) = length(g1.band) == 8192;

  g2 = g0;
  g2 = refine_cpgrid(g2, 1);
  g2 = refine_cpgrid(g2, 1);
  g2 = refine_cpgrid(g2, 1);
  g2 = refine_cpgrid(g2, 3);

  c = c + 1;
  pass(c) = length(g2.band) == length(g1.band);

