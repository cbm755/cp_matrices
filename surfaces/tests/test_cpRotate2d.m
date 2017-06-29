function [pass, str] = test_cpRotate2d()
  str = ['cpRotate2d test: rotate a line segment'];

  pass = [];
  c = 0;

  cp1 = @(x,y) cpLineSegment2d(x,y, 1/sqrt(2)*[-1.2 -1.2], 1/sqrt(2)*[1 1]);

  cpf = @(x,y) cpLineSegment2d(x,y, [-1.2 0], [1 0]);
  cp2 = @(x,y) cpRotate2d(x, y, cpf, pi/4);

  dx = 0.1;
  dy = 0.2;
  x1 = -2:dx:2;
  y1 = -2:dy:2;

  [x,y] = meshgrid(x1,y1);

  [cpx1 cpy1 dist1 bdy1] = cp1(x,y);
  [cpx2 cpy2 dist2 bdy2] = cp2(x,y);

  c=c+1;  pass(c) = assertAlmostEqual(cpx1, cpx2);
  c=c+1;  pass(c) = assertAlmostEqual(cpy1, cpy2);
  c=c+1;  pass(c) = assertAlmostEqual(dist1, dist2);
  c=c+1;  pass(c) = assertAlmostEqual(bdy1, bdy2);

  if ~(all(pass))
    max(max(abs(cpx1-cpx2)))
    max(max(abs(cpy1-cpy2)))
    max(max(abs(dist1-dist2)))
    max(max(abs(bdy1-bdy2)))

    figure(1);
    porcupine_plot2d(x, y, cpx1, cpy1, bdy1)
    figure(2)
    porcupine_plot2d(x, y, cpx2, cpy2, bdy2)
  end
