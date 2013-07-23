function [pass, str] = test_refine_grid3d_fixpt()
  str = ['refinement: fixed point iter converges to full grid, 3D'];
  pass = [];
  c = 0;

  g.dim = 3;
  g.dx = 0.25;
  g.x1d = (-2.0:g.dx:2.0)';
  g.y1d = g.x1d;
  g.z1d = g.x1d;
  [g.x, g.y, g.z] = meshgrid(g.x1d, g.y1d, g.z1d);
  g.x = g.x(:);
  g.y = g.y(:);
  g.z = g.z(:);
  g.cpfun = @cpSphere;
  g.band = 1:length(g.x);

  [g.cpx,g.cpy,g.cpz,g.dist] = g.cpfun(g.x, g.y, g.z);

  g1 = refine_cpgrid(g, 3);

  gA = g;
  % deliberately only give one CP
  gA.x = 1;
  gA.y = 0.5;
  gA.z = -0.3;
  [gA.cpx,gA.cpy,gA.cpz,gA.dist] = gA.cpfun(gA.x, gA.y, gA.z);

  g2 = refine_cpgrid(gA, 3);

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

  c = c + 1;
  pass(c) = length(g1.band) == 4784;

