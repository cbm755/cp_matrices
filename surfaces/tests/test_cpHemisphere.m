function [pass, str] = test_cpHemisphere()
  str = 'cpHemisphere tests, test a bunch of points';

  pass = [];
  c = 0;

  % right of origin
  [cpx,cpy,cpz,dist] = cpHemisphere(1,0,0);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy,cpz,dist], [1,0,0,0]);

  % below
  [cpx,cpy,cpz,dist] = cpHemisphere(1,0,-10);
  c = c + 1;
  pass(c) = assertAlmostEqual([cpx,cpy,cpz], [1,0,0]);

  % origin itself: don't care where it goes as long as its on the
  % surface
  [cpx,cpy,cpz,dist] = cpHemisphere(0,0,0);
  c = c + 1;
  pass(c) = assertAlmostEqual(norm([cpx,cpy,cpz],2), 1);
  c = c + 1;
  pass(c) = assertAlmostEqual(dist, 1);

  % origin again with radius 10
  [cpx,cpy,cpz,dist] = cpHemisphere(0,0,0,10);
  c = c + 1;
  pass(c) = assertAlmostEqual(norm([cpx,cpy,cpz],2), 10);
  c = c + 1;
  pass(c) = assertAlmostEqual(dist, 10);

