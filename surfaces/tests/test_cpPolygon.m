function [pass, str] = test_cpPolygon()
  str = ['cpPolygon tests: test a bunch of points, also sign of distance'];

  pass = [];
  c = 0;

  [cpx,cpy,sd] = cpPolygon(0,0);
  cp_tar = [0.05 0.2];
  sd_tar = -0.206155281280883;

  c = c + 1;
  pass(c) = sd < 0;

  c = c + 1;
  pass(c) = assertAlmostEqual([cpx cpy], cp_tar);
  c = c + 1;
  pass(c) = assertAlmostEqual(sd, sd_tar);


  [cpx,cpy,sd] = cpPolygon(10,10);
  c = c + 1;
  pass(c) = sd > 0;
  c = c + 1;
  pass(c) = sd > 10;

