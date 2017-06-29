function [pass, str] = test_cpRotate2d_arc()
  str = ['cpRotate2d test: rotate an arc'];

  pass = [];
  c = 0;

  rad = 0.95;

  for ii=1:4
    if ii==1, th = pi/10; end
    if ii==2, th = -pi/5; end
    if ii==3, th = 4.1*pi; end
    if ii==4, th = 0; end

    cp1 = @(x,y) cpArc(x,y, rad, [], th, th+pi);
    param1 = @(n) paramArc(n, rad, [], th, th+pi);

    cpf = @(x,y) cpSemicircle(x, y, rad);
    paramf = @(n) paramSemicircle(n, rad);

    cp2 = @(x,y) cpRotate2d(x, y, cpf, th);
    param2 = @(n) paramRotate2d(n, paramf, th);

    dx = 0.1;
    dy = 0.2;
    x1 = -2:dx:2;
    x1=setdiff(x1,0);  % remove origin
    y1 = -2:dy:2;
    [x,y] = meshgrid(x1,y1);

    [cpx1 cpy1 dist1 bdy1] = cp1(x,y);
    [cpx2 cpy2 dist2 bdy2] = cp2(x,y);

    c=c+1;  pass(c) = assertAlmostEqual(cpx1, cpx2);
    c=c+1;  pass(c) = assertAlmostEqual(cpy1, cpy2);
    c=c+1;  pass(c) = assertAlmostEqual(dist1, dist2);
    c=c+1;  pass(c) = assertAlmostEqual(bdy1, bdy2);

    if (1==0) ||  ~(all(pass))
      max(max(abs(cpx1-cpx2)))
      max(max(abs(cpy1-cpy2)))
      max(max(abs(dist1-dist2)))
      max(max(abs(bdy1-bdy2)))

      figure(1);
      porcupine_plot2d(x, y, cpx1, cpy1, bdy1)
      [xp,yp] = param1(256);  plot(xp,yp,'g-');
      figure(2)
      porcupine_plot2d(x, y, cpx2, cpy2, bdy2)
      [xp,yp] = param2(256);  plot(xp,yp,'g-');
      title(sprintf('cpRotate, \\theta=%g', th));
      pause
    end
  end

