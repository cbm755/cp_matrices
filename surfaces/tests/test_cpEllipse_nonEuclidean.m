function [pass, str] = test_cpEllipse_nonEuclidean()
  str = 'cpEllipse_nonEuclidean: test some points';

  pass = [];
  c = 0;

  a = 1.6;    b = 0.2;

  % check some points along the axes
  [cpx,cpy] = cpEllipse_nonEuclidean(2, 0, a, b);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [a 0], eps);

  [cpx,cpy] = cpEllipse_nonEuclidean(-2, 0, a, b);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [-a 0], eps);

  [cpx,cpy] = cpEllipse_nonEuclidean(0, 1, a, b);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [0 b], eps);

  [cpx,cpy] = cpEllipse_nonEuclidean(0, -1, a, b);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [0 -b], eps);

  % record the solution for a particular point
  [cpx,cpy] = cpEllipse_nonEuclidean(1, 1, a, b);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [0.971666837132019  0.15889577743848], 10*eps);

  % and should be exactly equal if we switch signs
  [cpx2,cpy2] = cpEllipse_nonEuclidean(-1, -1, a, b);
  c = c + 1;
  pass(c) = all([cpx cpy] == -[cpx2 cpy2]);

  % try a circle
  [cpx,cpy,sdist] = cpEllipse_nonEuclidean(1, 1, 1, 1);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [1 1]*(sqrt(2)/2), eps);
  c = c + 1;
  pass(c) = assertAlmostEqual(sdist, sqrt(2)-1, eps);

  [cpx,cpy,sdist] = cpEllipse_nonEuclidean(0.2, 0.2, 1, 1);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [1 1]*(sqrt(2)/2), eps);
  c = c + 1;
  pass(c) = assertAlmostEqual(sdist, sqrt(2*.2^2)-1, eps);

  % degenerate: lines
  [cpx,cpy] = cpEllipse_nonEuclidean(1, 1, 2, 0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [1 0], eps);

  [cpx,cpy] = cpEllipse_nonEuclidean(-1, 3, 2, 0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [-1 0], eps);

  [cpx,cpy] = cpEllipse_nonEuclidean(1, 0.1, 0, 1);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [0 0.1], eps);

