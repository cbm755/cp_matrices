function [pass, str] = test_weno6_2d_ellipse()
  str = 'test weno6 interpolation routines on an ellipse';

  c = 0;
  pass = [];

  dx = 0.14;
  dy = 0.1;
  maxdx = max([dx dy]);

  % make vectors of x, y, positions of the grid
  x1d = (-1.5-10*dx:dx:1.5+10*dx)';
  y1d = (-1-11*dy:dy:1+11*dy)';

  nx = length(x1d);
  ny = length(y1d);

  [xx yy] = meshgrid(x1d, y1d);
  [cpx, cpy, dist] = cpEllipse(xx,yy);
  cpxg = cpx(:); cpyg = cpy(:);


  %% Banding
  dim = 2;  % dimension
  p = 5;    % interpolation order
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(abs(dist) <= bw*maxdx);
  % store closest points in the band;
  cpxg = cpxg(band); cpyg = cpyg(band);
  xg = xx(band); yg = yy(band);

  % a smooth function on the comp domain (but not a surface function so
  % can't easily do accuracy test)
  u = cos(3*xg+2*yg.^3);

  %figure(1); clf; hold on;
  %plot2d_compdomain(u, xg, yg, dx, dy, 1)

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


  %% weno is "close" to quintic poly interp
  E = interp2_matrix(x1d, y1d, cpxg, cpyg, 5, band);
  wpoly = E*u;
  %max(abs(w1 - wpoly))
  %maxdx^6
  const = max(abs(w1-wpoly)) / maxdx^6;
  c = c + 1;
  pass(c) = (const < 200) & (const > .1);


  %% force ideal weights
  opt = [];
  opt.forceWeights = 10;
  opt.wenoeps = 1e-6;
  w3 = weno6_interp(cp, u, [cp.cpx cp.cpy], false, opt);


  c = c + 1;
  pass(c) = max(abs(w3 - wpoly)) < 25*eps;


  %% force ideal weights, middle stencil
  opt.forceWeights = 0;
  w3 = weno6_interp(cp, u, [cp.cpx cp.cpy], false, opt);
  E = interp2_matrix(x1d, y1d, cpxg, cpyg, 3, band);
  wpoly = E*u;
  %max(abs(w3 - wpoly))

  c = c + 1;
  pass(c) = max(abs(w3 - wpoly)) < 25*eps;

