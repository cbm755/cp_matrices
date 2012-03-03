function [pass, str] = test_weno6_2d()
  str = 'test weno6 interpolation routines 2D';

  c = 0;
  pass = [];

  dx = 0.08;
  dy = 0.05;
  maxdx = max([dx dy]);

  % make vectors of x, y, positions of the grid
  x1d = (-1-10*dx:dx:1+10*dx)';
  y1d = (-1-11*dy:dy:1+11*dy)';

  nx = length(x1d);
  ny = length(y1d);

  [xx yy] = meshgrid(x1d, y1d);
  [cpx, cpy, dist] = cpCircle(xx,yy);
  cpxg = cpx(:); cpyg = cpy(:);


  %% Banding
  dim = 2;  % dimension
  p = 5;    % interpolation order
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(abs(dist) <= bw*maxdx);
  % store closest points in the band;
  cpxg = cpxg(band); cpyg = cpyg(band);
  xg = xx(band); yg = yy(band);

  % a smooth function on the surface
  [th, r] = cart2pol(xg,yg);
  u = (cos(th)).^2;

  % grid object
  cp.x1d = x1d;
  cp.y1d = y1d;
  cp.cpx = cpxg;
  cp.cpy = cpyg;
  cp.dim = 2;
  cp.band = band;

  %% weno6 interp with and without caching
  w1 = weno6_interp(cp, u, [cp.cpx cp.cpy]);

  Weno6Cache = weno6_interp(cp, u, [cp.cpx cp.cpy], true);
  w2 = weno6_interp(Weno6Cache, u);

  c = c + 1;
  pass(c) = max(abs(w1-w2)) == 0;


  %% weno interp should be O(dx^{p+1}), p=5
  %max(abs(w1-u))
  %maxdx^6
  const = max(abs(w1-u)) / maxdx^6;
  c = c + 1;
  pass(c) = (const > 0.05) & (const < 20);
  %keyboard

  %% make sure the opt structure works
  opt = [];
  opt.forceWeights = -10;
  opt.wenoeps = 1e-6;
  w3 = weno6_interp(cp, u, [cp.cpx cp.cpy], false, opt);

  c = c + 1;
  pass(c) = max(abs(w1-w3)) == 0;


  %% force ideal weights
  opt.forceWeights = 10;
  w3 = weno6_interp(cp, u, [cp.cpx cp.cpy], false, opt);

  E = interp2_matrix(x1d, y1d, cpxg, cpyg, 5, band);
  wpoly = E*u;
  %max(abs(w3 - wpoly))

  c = c + 1;
  pass(c) = max(abs(w3 - wpoly)) < 25*eps;


  %% force ideal weights, middle stencil
  opt.forceWeights = 0;
  w3 = weno6_interp(cp, u, [cp.cpx cp.cpy], false, opt);
  E = interp2_matrix(x1d, y1d, cpxg, cpyg, 3, band);
  wpoly = E*u;
  %max(abs(w3 - wpoly))

  c = c + 1;
  pass(c) = max(abs(w3 - wpoly)) < 5*eps;

