function [pass, str] = test_cpPolygon_cpSquare_compare()
  str = ['cpPolygon test: compare to cpSquare'];

  pass = [];
  c = 0;

  poly = [-1    -1
           1    -1
           1     1
          -1     1];

  % delete the origin, cp is ambiguous there
  dx = 1;
  xx = [-3:dx:-dx  dx:dx:3];
  yy = [-2:dx:-dx  dx:dx:2];
  [x,y] = meshgrid(xx, yy);

  [cpx1, cpy1, sd1] = cpPolygon(x, y, poly);
  [cpx2, cpy2, sd2] = cpSquare(x, y);

  c = c + 1;
  pass(c) = assertAlmostEqual(cpx1, cpx2);
  c = c + 1;
  pass(c) = assertAlmostEqual(cpy1, cpy2);

  warning('two known failures: bug in cpPolygon');
  c = c + 1;
  pass(c) = assertAlmostEqual(sd1, sd2);


  % TODO: BUG: something wrong with the signed distance field in cpPolygon
  [cpx1, cpy1, sd1] = cpPolygon(-2, -1, poly)
  [cpx2, cpy2, sd2] = cpSquare(-2, -1)
  c = c + 1;
  pass(c) = assertAlmostEqual(sd1, 1);
  c = c + 1;
  pass(c) = assertAlmostEqual(sd2, 1);

  %pass
  %keyboard