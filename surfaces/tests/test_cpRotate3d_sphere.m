function [pass, str] = test_cpRotate3d_sphere()
  str = ['cpRotate3d test: rotate a sphere (no effect)'];

  c = 0;  pass = [];

  % random sets of three angles
  ths = pi*[2.6     2.3     0.1;  ...
           -0.6     3.6     3.1;  ...
            3.1    -1.4     0.7;  ...
           -0.9     1.4    -2.8;  ...
            2.2    -0.5    -2.4;  ...
           -0.8     2.7    -0.7;  ...
            2.5     2.2       2;  ...
           -2.3     3.9    -1.5];

  for i=1:size(ths,1)
    th = ths(i,:);

    cp1 = @(x,y,z) cpSphere(x,y,z);
    cp2 = @(x,y,z) cpRotate(x,y,z, cp1, th);

    dx = 0.2;  dy = 0.25;  dz = 0.5;
    x1 = -2:dx:2;  x1 = setdiff(x1,0);  % remove origin, ambiguous
    y1 = -2:dy:2;
    z1 = -2:dz:2;
    [x,y,z] = meshgrid(x1,y1,z1);

    [cpx1 cpy1 cpz1 dist1] = cp1(x,y,z);
    [cpx2 cpy2 cpz2 dist2] = cp2(x,y,z);

    c=c+1;  pass(c) = assertAlmostEqual(cpx1, cpx2) && ...
                      assertAlmostEqual(cpy1, cpy2) && ...
                      assertAlmostEqual(cpz1, cpz2) && ...
                      assertAlmostEqual(dist1, dist2);

    if ~(all(pass))
      max(max(max(abs(cpx1-cpx2))))
      max(max(max(abs(cpy1-cpy2))))
      max(max(max(abs(cpz1-cpz2))))
      max(max(max(abs(dist1-dist2))))
    end
  end
