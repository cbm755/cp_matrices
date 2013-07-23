function [pass, str] = test_refine_grid2d_fixpt()
  str = ['refinement: fixed point iter converges to full grid, 2D'];
  pass = [];
  c = 0;

  g.dim = 2;
  g.dx = 0.125;
  g.x1d = (-2.0:g.dx:2.0)';
  g.y1d = g.x1d;
  [g.x, g.y] = meshgrid(g.x1d, g.y1d);
  g.x = g.x(:);
  g.y = g.y(:);
  g.cpfun = @cpCircle;
  g.band = 1:length(g.x);


  [g.cpx,g.cpy,g.dist] = g.cpfun(g.x, g.y);

  g1 = refine_cpgrid(g, 3);

  gA = g;
  % deliberately only give one CP
  gA.x = 1;
  gA.y = 0;
  [gA.cpx,gA.cpy,gA.dist] = gA.cpfun(gA.x, gA.y);

  g2 = refine_cpgrid(gA, 3);

  %plot2d_compdomain(ones(size(g2.x)), g2.x, g2.y, g2.dx, g2.dx);

  c = c + 1;
  pass(c) = length(g1.band) == length(g2.band);
  c = c + 1;
  pass(c) = all(g1.dist == g2.dist);
  c = c + 1;
  pass(c) = all(g1.x == g2.x);
  c = c + 1;
  pass(c) = all(g1.y == g2.y);
  c = c + 1;
  pass(c) = all(g1.cpx == g2.cpx);
  c = c + 1;
  pass(c) = all(g1.cpy == g2.cpy);
