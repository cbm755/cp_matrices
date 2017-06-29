function [pass, str] = test_cpEllipseCurve()
  str = 'test some points for the (hacky) ellipseCurve';

  pass = [];
  c = 0;
  tol = 3*eps;

  % reasonable defaults
  [cpx,cpy] = cpEllipseCurve(0, 0);
  c = c + 1;
  pass(c) = assertAlmostEqual(cpx, 0, tol);


  a = 1.6;
  b = 0.2;
  lim = [pi/12  3*pi/4];
  leftx = a*cos(lim(2));
  lefty = b*sin(lim(2));
  rightx = a*cos(lim(1));
  righty = b*sin(lim(1));

  % pts on vertical axis
  [cpx,cpy] = cpEllipseCurve(0, 2, lim, [a b]);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [0 b], tol);

  [cpx,cpy] = cpEllipseCurve(0, 0, lim, [a b]);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [0 b], tol);

  % points to the right
  [cpx,cpy,dist,bdy] = cpEllipseCurve(2, 0, lim, [a b]);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [rightx righty], tol);
  c = c + 1;
  pass(c) = bdy == 2;
  [cpx,cpy,dist,bdy] = cpEllipseCurve(2, 1, lim, [a b]);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [rightx righty], tol);
  c = c + 1;
  pass(c) = bdy == 2;

  % points to the left
  [cpx,cpy,dist,bdy] = cpEllipseCurve(-3, 0, lim, [a b]);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [leftx lefty], tol);
  c = c + 1;
  pass(c) = bdy == 1;

  [cpx,cpy,dist,bdy] = cpEllipseCurve(-2.5, 5, lim, [a b]);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], [leftx lefty], tol);
  c = c + 1;
  pass(c) = bdy == 1;

  % points on the x axis
  x = 6*rand(10,1) - 4;
  y = zeros(size(x));
  [cpx,cpy] = cpEllipseCurve(x, y, lim, [a b]);
  c = c + 1;
  pass(c) = all(cpy > 0);

  % test some random points to see if dist matches
  x = 6*rand(100,1) - 3;
  y = 3*rand(100,1);
  [cpx,cpy,dist] = cpEllipseCurve(x, y, lim, [a b]);
  c = c + 1;
  dist2 = sqrt((cpx - x).^2 + (cpy - y).^2);
  pass(c) = assertAlmostEqual(dist,dist2);
  %I = abs(dist - dist2) > eps;
  %figure(1); clf;
  %porcupine_plot2d(x,y,cpx,cpy)
  %plot(x(I), y(I), 'ro');

  % compare to circle
  x = 6*rand(100,1) - 3;
  y = 3*rand(100,1);
  [cpx,cpy,dist] = cpEllipseCurve(x, y, [0 pi], [1 1]);
  [cpx2,cpy2,dist2] = cpSemicircle(x, y, 1);
  c = c + 1;
  pass(c) = assertAlmostEqual(cpx,cpx2);
  c = c + 1;
  pass(c) = assertAlmostEqual(cpy,cpy2);
  c = c + 1;
  pass(c) = assertAlmostEqual(dist,dist2);

