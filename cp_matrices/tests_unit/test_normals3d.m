function [pass, str] = test_normals3d()
  str = 'normals test 3d: recover normals from cp-rep';


  %% Construct a grid in the embedding space
  dx = 0.2;
  x1d = (-2:dx:2)';
  y1d = x1d;
  z1d = x1d;


  %% Find closest points on the surface
  [xx yy zz] = meshgrid(x1d, y1d, z1d);
  [cpx cpy cpz sdist] = cpSphere(xx,yy,zz);
  %[cpx cpy cpz sdist] = cpEllipsoid(xx,yy,zz, [1.3 0.4]);


  %% Banding
  dim = 3;
  p = 3;
  order = 2;
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
  band = find(abs(sdist) <= bw*dx);
  % store closest points in the band;
  cpxg = cpx(band); cpyg = cpy(band); cpzg = cpz(band);
  xg = xx(band); yg = yy(band); zg = zz(band);
  sdistg = sdist(band);


  %% tests
  [n1,n2,n3] = normals_from_cp(xx,yy,zz, cpx,cpy,cpz, ...
                               sdist, dx, 0);

  % where is origin?
  O = (xx == 0 & yy == 0 & zz == 0);

  % origin should have reasonable normal
  test1 = assertAlmostEqual(n1(O).^2 + n2(O).^2 + n3(O).^2, 1);
  pass = [test1];

  % other than origin, should give close to exact normal
  A = xx ./ sqrt(xx.^2 + yy.^2 + zz.^2);
  B = yy ./ sqrt(xx.^2 + yy.^2 + zz.^2);
  C = zz ./ sqrt(xx.^2 + yy.^2 + zz.^2);
  test1 = assertAlmostEqual(n1(~O), A(~O), 1000*eps);
  test2 = assertAlmostEqual(n2(~O), B(~O), 1000*eps);
  test3 = assertAlmostEqual(n3(~O), C(~O), 1000*eps);
  pass = [pass test1 test2 test3];


  %% now test with bands
  [n1g,n2g,n3g] = normals_from_cp(xg,yg,zg, cpxg,cpyg,cpzg, ...
                                  sdistg, dx, 0);

  A = xg ./ sqrt(xg.^2 + yg.^2 + zg.^2);
  B = yg ./ sqrt(xg.^2 + yg.^2 + zg.^2);
  C = zg ./ sqrt(xg.^2 + yg.^2 + zg.^2);
  test1 = assertAlmostEqual(n1g, A, 1e-2);
  test2 = assertAlmostEqual(n2g, B, 1e-2);
  test3 = assertAlmostEqual(n3g, C, 1e-2);
  pass = [pass test1 test2 test3];

  % the ones that require such a big tolerance are (0.8,0.6,0) and
  % similar
  %I = find(abs(n1g-A) > 1e-10);
  %[xg(I) yg(I) zg(I)]

  %% now without orientation
  [n1g,n2g,n3g] = normals_from_cp(xg,yg,zg, cpxg,cpyg,cpzg, ...
                                  abs(sdistg), dx, 0);

  A = xg ./ sqrt(xg.^2 + yg.^2 + zg.^2);
  B = yg ./ sqrt(xg.^2 + yg.^2 + zg.^2);
  C = zg ./ sqrt(xg.^2 + yg.^2 + zg.^2);
  test1 = all((n1g - A) < 1e-2 | (n1g + A) < 1e-2);
  test2 = all((n2g - B) < 1e-2 | (n2g + B) < 1e-2);
  test3 = all((n3g - C) < 1e-2 | (n3g + C) < 1e-2);
  pass = [pass test1 test2 test3];


