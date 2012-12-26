function [pass, str] = test_orientation1()
  str = 'orientation test: recover orientation from cp-rep and seeds';


  %% Construct a grid in the embedding space
  dx = 0.2;
  x1d = (-2:dx:2)';
  y1d = x1d;
  z1d = x1d;


  %% Find closest points on the surface
  [xx yy zz] = meshgrid(x1d, y1d, z1d);
  %[cpx cpy cpz sdist] = cpSphere(xx,yy,zz);
  [cpx cpy cpz sdist] = cpEllipsoid(xx,yy,zz, [1.3 0.4]);
  % make into vectors
  cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);


  %% Banding
  dim = 3;
  p = 3;
  order = 2;
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
  band = find(abs(sdist) <= bw*dx);
  % store closest points in the band;
  cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
  xg = xx(band); yg = yy(band); zg = zz(band);
  sdistg = sdist(band);


  %% we do some tests on the individual stages
  seedin = find((abs(xx - 2*dx) < 10*eps) & yy == 0 & zz == 0);
  test1 = sdist(seedin) < 0;
  inside = orientation_fill(xx,yy,zz,sdist,dx,seedin,0);
  test2 = all(sdist(logical(inside)) < 0);

  seedout = find((abs(xx - x1d(end)) < 10*eps) & yy == 0 & zz == 0);
  test3 = sdist(seedout) > 0;
  outside = orientation_fill(xx,yy,zz,sdist,dx,seedout,0);
  test4 = all(sdist(logical(outside)) > 0);

  pass = [test1 test2 test3 test4];

  % now some counting checks
  ni = nnz(inside);
  no = nnz(outside);
  total = length(xx(:));
  test1 = (ni + no <= total);

  noclass = find(~inside & ~outside);
  nnc = length(noclass);
  test2 = ni + no + length(noclass) == total;
  pass = [pass test1 test2];

  % check if nonclassified are at least close to surface
  % TODO: probably will fail until we tweak parameters in _fill
  d = abs(sdist(noclass));
  test = all(d < dx);
  pass = [pass test];


  %% stage 2 processing
  [in2,out2] = orientation_stage2(xx,yy,zz,cpx,cpy,cpz,sdist,dx,inside,outside,0);

  ni2 = nnz(in2);
  no2 = nnz(out2);
  total = length(xx(:));

  noclass2 = find(~in2 & ~out2);
  nnc2 = length(noclass2);
  test1 = nnc2 == 0;
  test2 = ni2 + no2 == total;

  pass = [pass test1 test2];


  %% now use the main interface
  [in3,sdist2] = orientation_from_cp(xx,yy,zz, cpx,cpy,cpz, ...
                                     sdist, dx, ...
                                     seedin, seedout);

  test1 = all(all(all(in2 == in3)));
  test2 = all(all(all(sdist == sdist2)));  % warning FP equality test

  pass = [pass test1 test2];


  %% now test with bands
  seedin = find((abs(xg - 2*dx) < 10*eps) & yg == 0 & zg == 0);
  seedout = find((abs(xg - x1d(end)) < 10*eps) & yg == 0 & zg == 0);

  [ing,sdistg2] = orientation_from_cp(xg,yg,zg, cpxg,cpyg,cpzg, ...
                                      sdistg, dx, ...
                                      seedin, seedout, 0);

  test = all(sdistg == sdistg2);  % warning FP equality test

  pass = [pass test];

