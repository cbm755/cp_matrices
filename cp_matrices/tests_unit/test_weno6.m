function [pass, str] = test_weno6()
  str = 'test weno6 interpolation routines 3D';

  c = 0;
  pass = [];

  dx = 0.2;
  dy = 0.25;
  dz = 0.3;
  maxdx = max([dx dy dz]);

  % make vectors of x, y, positions of the grid
  x1d = (-1-10*dx:dx:1+10*dx)';
  y1d = (-1-11*dy:dy:1+11*dy)';
  z1d = (-1-12*dz:dz:1+12*dz)';

  nx = length(x1d);
  ny = length(y1d);
  nz = length(z1d);

  [xx yy zz] = meshgrid(x1d, y1d, z1d);
  [cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
  cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);


  %% Banding
  dim = 3;  % dimension
  p = 5;    % interpolation order
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(abs(dist) <= bw*maxdx);
  % store closest points in the band;
  cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
  xg = xx(band); yg = yy(band); zg = zz(band);

  % a smooth function on the surface
  [th, phi, r] = cart2sph(xg,yg,zg);
  u = (cos(phi + pi/2)).^2;

  % grid object
  cp.x1d = x1d;
  cp.y1d = y1d;
  cp.z1d = z1d;
  cp.cpx = cpxg;
  cp.cpy = cpyg;
  cp.cpz = cpzg;
  cp.dim = 3;
  cp.band = band;

  %% weno6 interp with and without caching
  w1 = weno6_interp(cp, u, [cp.cpx cp.cpy cp.cpz]);

  Weno6Cache = weno6_interp(cp, u, [cp.cpx cp.cpy cp.cpz], true);
  w2 = weno6_interp(Weno6Cache, u);

  c = c + 1;
  pass(c) = max(abs(w1-w2)) == 0;


  %% weno interp should be O(dx^{p+1}), p=5
  %max(abs(w1-u))
  %maxdx^6
  const = max(abs(w1-u)) / maxdx^6;
  c = c + 1;
  pass(c) = (const > 0.05) & (const < 20);


  %% make sure the opt structure works
  opt = [];
  opt.forceWeights = -10;
  opt.wenoeps = 1e-6;
  w3 = weno6_interp(cp, u, [cp.cpx cp.cpy cp.cpz], false, opt);

  c = c + 1;
  pass(c) = max(abs(w1-w3)) == 0;


  %% force ideal weights
  opt.forceWeights = 10;
  w3 = weno6_interp(cp, u, [cp.cpx cp.cpy cp.cpz], false, opt);

  E = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, 5, band);
  wpoly = E*u;
  %max(abs(w3 - wpoly))

  c = c + 1;
  pass(c) = max(abs(w3 - wpoly)) < 25*eps;


  %% force ideal weights, middle stencil
  opt.forceWeights = 0;
  w3 = weno6_interp(cp, u, [cp.cpx cp.cpy cp.cpz], false, opt);
  E = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, 3, band);
  wpoly = E*u;
  %max(abs(w3 - wpoly))

  c = c + 1;
  pass(c) = max(abs(w3 - wpoly)) < 5*eps;

