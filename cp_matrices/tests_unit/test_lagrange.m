function [pass, str] = test_lagrange()
  str = 'some fixed cases to test the Barycentric Lagrange weights';

  c = 0;
  pass = [];

  % first 4 are just binary weights testing
  xg = 0;  x = 0;  dx = 1;  N = 4;  w = [1 0 0 0];
  w1 = LagrangeWeights1D_vec(xg,x,dx,N);
  c = c + 1;
  pass(c) = all(w1 == w);

  xg = 0;  x = 1;  dx = 1;  N = 4;  w = [0 1 0 0];
  w1 = LagrangeWeights1D_vec(xg,x,dx,N);
  c = c + 1;
  pass(c) = all(w1 == w);

  xg = 0;  x = 2;  dx = 1;  N = 4;  w = [0 0 1 0];
  w1 = LagrangeWeights1D_vec(xg,x,dx,N);
  c = c + 1;
  pass(c) = all(w1 == w);

  xg = 0;  x = 3;  dx = 1;  N = 4;  w = [0 0 0 1];
  w1 = LagrangeWeights1D_vec(xg,x,dx,N);
  c = c + 1;
  pass(c) = all(w1 == w);

  % some special cases
  xg = 0;  x = 1.5;  dx = 1;  N = 4;  w = [-1 9 9 -1]/16;
  w1 = LagrangeWeights1D_vec(xg,x,dx,N);
  c = c + 1;
  pass(c) = all(w1 == w);

  xg = 0;  x = 1.25;  dx = 1;  N = 4;  w = [-7  105  35  -5]/128;
  w1 = LagrangeWeights1D_vec(xg,x,dx,N);
  c = c + 1;
  pass(c) = assertAlmostEqual(w1, w);

  % nearest neighbour
  xg = 0;  x = 0.5;  dx = 1;  N = 1;  w = 1;
  w1 = LagrangeWeights1D_vec(xg,x,dx,N);
  c = c + 1;
  pass(c) = assertAlmostEqual(w1, w);

