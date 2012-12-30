function [pass, str] = test_normals2d()
  str = 'normals test 2d: recover normals from cp-rep';


  %% Construct a grid in the embedding space
  dx = 0.1;
  x1d = (-2:dx:2)';
  y1d = x1d;


  %% Find closest points on the surface
  [xx yy] = meshgrid(x1d, y1d);
  [cpx cpy sdist] = cpCircle(xx, yy);
  paramf = @paramCircle;
  %[cpx cpy sdist] = cpEllipse(xx,yy, 1.3, 0.2);
  %paramf = @(n) paramEllipse(n, 1.3, 0.2);

  %% Banding
  dim = 2;
  p = 3;
  order = 2;
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
  band = find(abs(sdist) <= bw*dx);
  % store closest points in the band;
  cpxg = cpx(band); cpyg = cpy(band);
  xg = xx(band); yg = yy(band);
  sdistg = sdist(band);


  %% tests
  [n1,n2] = normals_from_cp(xx,yy, cpx,cpy, sdist, dx, 0);

  % where is origin?
  O = (xx == 0 & yy == 0);

  % origin should have reasonable normal
  test1 = assertAlmostEqual(n1(O).^2 + n2(O).^2, 1);
  pass = [test1];

  % other than origin, should give close to exact normal
  A = xx ./ sqrt(xx.^2 + yy.^2);
  B = yy ./ sqrt(xx.^2 + yy.^2);
  test1 = assertAlmostEqual(n1(~O), A(~O), 1000*eps);
  test2 = assertAlmostEqual(n2(~O), B(~O), 1000*eps);
  pass = [pass test1 test2];


  %% now test with bands, oriented
  [n1g,n2g] = normals_from_cp(xg,yg, cpxg,cpyg, sdistg, dx, 0);

  A = xg ./ sqrt(xg.^2 + yg.^2);
  B = yg ./ sqrt(xg.^2 + yg.^2);
  test1 = assertAlmostEqual(n1g, A, 1e-3);
  test2 = assertAlmostEqual(n2g, B, 1e-3);
  pass = [pass test1 test2];

  % the ones that need such a large tolerances are at (0.6,0.8)
  %I = find(abs(n1g-A) > 1e-10)
  %keyboard


  %% now test with bands, non-oriented
  [n1g,n2g] = normals_from_cp(xg,yg, cpxg,cpyg, abs(sdistg), dx, 0);

  A = xg ./ sqrt(xg.^2 + yg.^2);
  B = yg ./ sqrt(xg.^2 + yg.^2);
  test1 = all((n1g - A) < 1e-3 | (n1g + A) < 1e-3);
  test2 = all((n2g - B) < 1e-3 | (n2g + B) < 1e-3);
  pass = [pass test1 test2];

