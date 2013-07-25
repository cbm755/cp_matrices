function [pass, str] = test_monocubic()
  str = 'test monocubic interp in 2d (data bounded, others?)';

  c = 0;
  pass = [];

  dx = 0.1;
  dy = 0.15;
  maxdx = max([dx dy]);

  % make vectors of x, y, positions of the grid
  x1d = (-1-6*dx:dx:1+6*dx)';
  y1d = (-1-7*dy:dy:1+7*dy)';

  nx = length(x1d);
  ny = length(y1d);

  [xx yy] = meshgrid(x1d, y1d);
  [cpx, cpy, dist] = cpCircle(xx,yy);
  cpxg = cpx(:); cpyg = cpy(:);


  %% Banding
  dim = 2;  % dimension
  p = 3;    % interpolation degree (for MC, want 3 in below formula)
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(abs(dist) <= bw*maxdx);
  cpxg = cpxg(band); cpyg = cpyg(band);
  xg = xx(band); yg = yy(band);

  % some smooth function on the surface
  [th, r] = cart2pol(xg,yg);
  u = (cos(th + pi/2)).^2;

  % random values
  u_rand = rand(size(xg));

  % grid object
  cp.x1d = x1d;
  cp.y1d = y1d;
  cp.cpx = cpxg;
  cp.cpy = cpyg;
  cp.dim = 2;
  cp.band = band;

  % data bounded: random values test
  w1 = monocubic_interp(cp, u_rand, [cp.cpx cp.cpy]);
  c = c + 1;
  pass(c) = min(w1) >= min(u_rand);
  c = c + 1;
  pass(c) = max(w1) <= max(u_rand);

  % TODO: a timing test?
  %T1 = cputime();
  %w1 = monocubic_interp(cp, u, [cp.cpx cp.cpy]);
  %T1 = cputime()-T1;

  %E = interp2_matrix(x1d, y1d, cpxg, cpyg, 2, band);
  %w2 = E*u;
  %max(abs(w2-u))

