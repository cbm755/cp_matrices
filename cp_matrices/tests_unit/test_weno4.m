function [pass, str] = test_weno4()
  str = 'test weno4 interpolation routines';

  c = 0;
  pass = [];

  warning('TODO: fix the weno4 code to use different dx,dy,dz, relpt')
  dx = 0.1;
  dy = 0.1;
  dz = 0.1;
  maxdx = max([dx dy dz]);

  % make vectors of x, y, positions of the grid
  x1d = (-1-5*dx:dx:1+5*dx)';
  y1d = (-1-5*dy:dy:1+5*dy)';
  z1d = (-1-5*dz:dz:1+5*dz)';

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
  %T = tic;
  w1 = weno4_interp(cp, u, [cp.cpx cp.cpy cp.cpz]);
  %toc(T)

  disp('*** weno4_caching call but without caching');
  %T = tic;
  w2 = weno4_interp_caching(cp, u, [cp.cpx cp.cpy cp.cpz]);
  %toc(T)

  c = c + 1;
  pass(c) = max(abs(w1-w2)) == 0;


  %disp('*** weno4_caching: calling build cache');
  WenoCache = weno4_interp_caching(cp, u, [cp.cpx cp.cpy cp.cpz], 'cache');

  %disp('*** cache built, now interpolating');
  w3 = weno4_interp_caching(WenoCache, u);
  w4 = weno4_interp_caching(WenoCache, 1.1*u);
  w5 = weno4_interp_caching(WenoCache, 1.2*u);

  c = c + 1;
  pass(c) = max(abs(w1-w3)) == 0;

  %disp('should be O(dx^{p+1}), p=3')
  %max(abs(w1-u))
  %dx^4
  const = max(abs(w1-u)) / maxdx^4;
  c = c + 1;
  pass(c) = (const > 0.2) & (const < 5);

