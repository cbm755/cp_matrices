function [pass, str] = test_weno4()
  str = 'test weno4 interpolation routines (3D)';

  c = 0;
  pass = [];

  warning('TODO: fix the weno4 code to use different dx,dy,dz, relpt')
  dx = 0.1;
  dy = 0.12;
  dz = 0.14;
  maxdx = max([dx dy dz]);

  % make vectors of x, y, positions of the grid
  x1d = (-1-6*dx:dx:1+6*dx)';
  y1d = (-1-7*dy:dy:1+7*dy)';
  z1d = (-1-8*dz:dz:1+8*dz)';

  nx = length(x1d);
  ny = length(y1d);
  nz = length(z1d);

  [xx yy zz] = meshgrid(x1d, y1d, z1d);
  [cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
  cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);


  %% Banding
  dim = 3;  % dimension
  p = 3;    % interpolation order
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(abs(dist) <= bw*maxdx);
  % store closest points in the band;
  cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
  xg = xx(band); yg = yy(band); zg = zz(band);

  % some smooth function on the surface
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

  disp('*** weno4 call');
  % do it twice as first one might not give good timing
  w1 = weno4_interp(cp, u, [cp.cpx cp.cpy cp.cpz]);
  T1 = cputime();
  w1 = weno4_interp(cp, u, [cp.cpx cp.cpy cp.cpz]);
  T1 = cputime()-T1;

  disp('*** weno4_caching: calling build cache');
  WenoCache = weno4_interp(cp, u, [cp.cpx cp.cpy cp.cpz], 'cache');
  T2 = cputime();
  WenoCache = weno4_interp(cp, u, [cp.cpx cp.cpy cp.cpz], 'cache');
  T2 = cputime()-T2;

  disp('*** cache built, now interpolating');
  w3 = weno4_interp(WenoCache, u);
  T3 = cputime();
  w3 = weno4_interp(WenoCache, u);
  T3 = cputime()-T3;
  %w4 = weno4_interp(WenoCache, 1.1*u);
  %w5 = weno4_interp(WenoCache, 1.2*u);

  %% with or without caching gives same result
  c = c + 1;
  pass(c) = max(abs(w1-w3)) == 0;

  %% one call without caching should be faster
  % this is dangerous b/c non-deterministic
  %c = c + 1;
  %pass(c) = (T1/(T2+T3) < 0.75);


  %% subsequent calls should be faster because of caching
  % another non-deterministic
  %c = c + 1;
  %pass(c) = (T3/T1 < 0.8);


  c = c + 1;
  pass(c) = max(abs(w1-w3)) == 0;

  %disp('should be O(dx^{p+1}), p=3')
  %max(abs(w1-u))
  %dx^4
  const = max(abs(w1-u)) / maxdx^4;
  c = c + 1;
  pass(c) = (const > 0.2) & (const < 5);

