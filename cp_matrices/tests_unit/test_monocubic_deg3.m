function [pass, str] = test_monocubic_deg3()
  str = 'test monocubic and other interpolation degrees in 2D';

  c = 0;
  pass = [];

  for i=1:3
    [errMC(i) errQ(i) errC(i)] = helper(i);
  end

  % quadratic p=2 should be third-order
  A = errMC(1:end-1) ./ errMC(2:end);
  c = c + 1;
  pass(c) = all(A > 5);  % should be 8

  % quadratic p=2 should be third-order
  A = errQ(1:end-1) ./ errQ(2:end);
  c = c + 1;
  pass(c) = all(A > 7);  % should be 8

  % cubic p=3 should be fourth-order
  A = errC(1:end-1) ./ errC(2:end);
  c = c + 1;
  pass(c) = all(A > 14); % should be 16
end


function [err err2 err3] = helper(n)
  dx = 0.2/2^(n);
  dy = 0.15/2^(n);
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
  p = 3;    % interpolation order
  bw = rm_bandwidth(dim, p);
  band = find(abs(dist) <= bw*maxdx);
  % store closest points in the band;
  cpxg = cpxg(band); cpyg = cpyg(band);
  xg = xx(band); yg = yy(band);

  % some smooth function on the surface
  [th, r] = cart2pol(xg,yg);
  u = (cos(th + pi/2)).^2;
  %u = rand(size(u));

  % grid object
  cp.x1d = x1d;
  cp.y1d = y1d;
  cp.cpx = cpxg;
  cp.cpy = cpyg;
  cp.dim = 2;
  cp.band = band;

  % do it twice as first one might not give good timing
  %w1 = monocubic_interp(cp, u, [cp.cpx cp.cpy]);
  T1 = cputime();
  w1 = monocubic_interp(cp, u, {cp.cpx cp.cpy});
  T1 = cputime()-T1;

  err = max(abs(w1-u));

  E = interp2_matrix(x1d, y1d, cpxg, cpyg, 2, band);
  w2 = E*u;
  err2 = max(abs(w2-u));

  E = interp2_matrix(x1d, y1d, cpxg, cpyg, 3, band);
  w3 = E*u;
  err3 = max(abs(w3-u));

end
