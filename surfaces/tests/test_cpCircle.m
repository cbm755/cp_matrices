function [pass, str] = test_cpCircle()
  str = 'cpCircle tests, test a bunch of points';

  pass = [];
  c = 0;

  % right of origin
  [cpx,cpy,dist] = cpCircle(1,0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy,dist], [1,0,0]);

  [cpx,cpy] = cpCircle(0.1,0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy], [1,0]);
  [cpx,cpy] = cpCircle(2,0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy], [1,0]);

  % just left of origin
  [cpx,cpy] = cpCircle(-0.1,0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy], [-1,0]);
  [cpx,cpy] = cpCircle(-1000,0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy], [-1,0]);

  % origin itself: don't care where it goes as long as its on the
  % surface
  [cpx,cpy,dist] = cpCircle(0,0);
  c = c + 1;
  pass(c) = assertAlmostEqual(cpx^2+cpy^2, 1);
  c = c + 1;
  pass(c) = assertAlmostEqual(dist, -1);

  % origin again with radius 10
  [cpx,cpy,sdist] = cpCircle(0,0,10);
  c = c + 1;
  pass(c) = assertAlmostEqual(sqrt(cpx^2+cpy^2), 10);
  c = c + 1;
  pass(c) = assertAlmostEqual(sdist, -10);



  [cpx,cpy] = cpCircle(2,2);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy], [sqrt(2)/2,sqrt(2)/2]);

  [cpx,cpy] = cpCircle(2,2,3);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy], [3*sqrt(2)/2,3*sqrt(2)/2]);

  [cpx,cpy] = cpCircle(2,2,3,[2,2]);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy], [2+3, 2]);

  [cpx,cpy] = cpCircle(10,10,3,[2,2]);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy], [2+3*sqrt(2)/2,2+3*sqrt(2)/2]);
