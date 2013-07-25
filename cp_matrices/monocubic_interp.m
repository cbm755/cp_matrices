function w = monocubic_interp(cp, f, xyz)
%MONOCUBIC_INTERP
%
%

%  monocubic_interp_2d()
%function w = monocubic_interp_2d(cp, )
  tic
  x1d = cp.x1d;
  y1d = cp.y1d;
  Nx = length(x1d);
  Ny = length(y1d);
  dx = x1d(2) - x1d(1);
  dy = y1d(2) - y1d(1);

  relpt = [cp.x1d(1)  cp.y1d(1)];

  if iscell(xyz)
    x = xyz{1};
    y = xyz{2};
  else
    x = xyz(:,1);
    y = xyz(:,2);
  end

  % determine the basepoint, roughly speaking this is "floor(xy)"
  % in terms of the grid
  [ijk,X] = findGridInterpBasePt_vec({x y}, 3, relpt, [dx dy]);
  % +1 here because the basepoint is actually the lowerleft corner
  % of the stencil and we want the "floor".
  xi = X{1} + dx;
  yi = X{2} + dy;
  ijk{1} = ijk{1} + 1;
  ijk{2} = ijk{2} + 1;

  I = sub2ind([Ny Nx], ijk{2}, ijk{1});
  B = findInBand(I, cp.band, Nx*Ny);

  [E W N S] = neighbourMatrices(cp, cp.band, cp.band);
  preptime = toc;

  opt = [];
  tic
  g = S*f;     u1 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = f;       u2 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = N*f;     u3 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = N*(N*f); u4 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);

  w = helper1d(u1, u2, u3, u4, yi, dy, y, opt);
  %w(i) = pchip(xi(i)+[0 1 2 3]*dy, [u1(i) u2(i) u3(i) u4(i)], yi(i));
  wenotime = toc;
  fprintf('monocubic: preptime=%g, calctime=%g\n', preptime, wenotime)
end


function u = helper1d(fim1, fi, fip1, fip2, xi, dx, x, opt)
%HELPER1D

  % this would return a matrix: each interpolant evaluated at each CP.
  % Anyway, this doesn't work either.  One approach is to implement a
  % vectorized version from scr
  U = pchip([-dx  0  dx  2*dx], [fim1 fi fip1 fip2], x-xi);
  u = diag(U);
end

function u = helper1d_slow(fim1, fi, fip1, fip2, xi, dx, x, opt)
%HELPER1D

  u = zeros(size(x));
  for i=1:length(x)
    u(i) = pchip(xi(i)+[-1 0 1 2]*dx, [fim1(i) fi(i) fip1(i) fip2(i)], x(i));
  end
end
