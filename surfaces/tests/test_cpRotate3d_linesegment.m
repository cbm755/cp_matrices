function [pass, str] = test_cpRotate3d_linesegment()
  str = ['cpRotate3d test: rotate a line segment'];

  c = 0;  pass = [];

  for i=1:7
    if i==1, th = [0 0 0];       p = [-1.2 0 0]; q = [1 0 0]; end
    if i==2, th = [pi 0 0];      p = [1.2 0 0];  q = [-1 0 0]; end
    if i==3, th = [pi/2 0 0];    p = [0 -1.2 0]; q = [0 1 0]; end
    if i==4, th = [0 pi/2 0];    p = [0 0 -1.2]; q = [0 0 1]; end
    if i==5, th = [0 0 pi/2];    p = [-1.2 0 0]; q = [1 0 0]; end
    if i==6, th = [pi/2 0 pi/2]; p = [0 0 -1.2]; q = [0 0 1]; end
    if i==7
      th = [pi/4 0 pi/4];
      p = [-1.2/sqrt(2) -1.2/2 -1.2/2];
      q = [1/sqrt(2) 1/2 1/2];
    end

    cpf = @(x,y,z) cpLineSegment3d(x,y,z, [-1.2 0 0], [1 0 0]);
    paramf = @(n) paramLineSegment3d(n, [-1.2 0 0], [1 0 0]);
    cp1 = @(x,y,z) cpRotate(x,y,z, cpf, th);
    param1 = @(n) paramRotate(n, paramf, th);
    cp2 = @(x,y,z) cpLineSegment3d(x,y,z, p, q);
    param2 = @(n) paramLineSegment3d(n, p, q);

    dx = 0.25;  dy = 0.5;  dz = 1;
    x1 = -2:dx:2;
    y1 = -2:dy:2;
    z1 = -2:dz:2;
    [x,y,z] = meshgrid(x1,y1,z1);

    [cpx1 cpy1 cpz1 dist1 bdy1] = cp1(x,y,z);
    [cpx2 cpy2 cpz2 dist2 bdy2] = cp2(x,y,z);

    c=c+1;  pass(c) = assertAlmostEqual(cpx1, cpx2) && ...
                      assertAlmostEqual(cpy1, cpy2) && ...
                      assertAlmostEqual(cpz1, cpz2) && ...
                      assertAlmostEqual(dist1, dist2);
    c=c+1;  pass(c) = assertAlmostEqual(bdy1, bdy2);

    if (1==0) || ~(all(pass))
      max(max(max(abs(cpx1-cpx2))))
      max(max(max(abs(cpy1-cpy2))))
      max(max(max(abs(cpz1-cpz2))))
      max(max(max(abs(dist1-dist2))))
      max(max(max(abs(bdy1-bdy2))))

      figure(1);
      porcupine_plot3d(x,y,z, cpx1,cpy1,cpz1, bdy1, param1)
      title([ 'cpRotate, angles=' num2str(th) ]);
      figure(2);
      porcupine_plot3d(x,y,z, cpx2,cpy2,cpz2, bdy2, param2)
      pause
    end
  end
