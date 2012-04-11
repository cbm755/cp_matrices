function [pass, str] = test_fd_weno5_1d()
  str = 'tests of finite difference WENO procedure in 1D';

  pass = [];
  c = 0;

  v = [1 1 1 1 1];
  y = fd_weno5_1d(v);

  c = c + 1;
  pass(c) = assertAlmostEqual(y, 1);

  y = fd_weno5_1d(10*v);
  c = c + 1;
  pass(c) = assertAlmostEqual(y, 10);