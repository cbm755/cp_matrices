function [pass, str] = test_cpSurfOfRevolution()
  str = 'cpSurfOfRevolution, rotate a circle to make a sphere';

  pass = [];
  c = 0;

  % delete the origin, cp is ambiguous there
  dx = 0.1;
  xx = [-2:dx:-dx  dx:dx:2];
  yy = [-2:dx:-dx  dx:dx:2];
  zz = [-2:dx:-dx  dx:dx:2];
  [x,y,z] = meshgrid(xx, yy, zz);

  [cpx1 cpy1 cpz1 sd1] = cpSphere(x, y, z);
  [cpx2 cpy2 cpz2 sd2] = cpSurfOfRevolution(x, y, z, @cpCircle);

  % these are set to pass with zero error
  c = c + 1;
  pass(c) = assertAlmostEqual(cpx1, cpx2, 5*eps);
  c = c + 1;
  pass(c) = assertAlmostEqual(cpy1, cpy2, 5*eps);
  c = c + 1;
  pass(c) = assertAlmostEqual(cpz1, cpz2, 5*eps);
  c = c + 1;
  pass(c) = assertAlmostEqual(sd1, sd2, 5*eps);



  [cpx1 cpy1 cpz1 sd1] = cpSphere(0,0,0);
  [cpx2 cpy2 cpz2 sd2] = cpSurfOfRevolution(0, 0, 0, @cpCircle);

  % no particular reason these should pass at all, origin is ambiguous
  c = c + 1;
  pass(c) = assertAlmostEqual(cpx1, cpx2, 0);
  c = c + 1;
  pass(c) = assertAlmostEqual(cpy1, cpy2, 0);
  c = c + 1;
  pass(c) = assertAlmostEqual(cpz1, cpz2, 0);
  c = c + 1;
  pass(c) = assertAlmostEqual(sd1, sd2, 0);

