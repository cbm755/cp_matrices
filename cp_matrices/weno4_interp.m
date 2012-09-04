function w = weno4_interp(cp, f, x, makeCache, opt)
%WENO4_INTERP  nonlinear WENO interpolation in 2D/3D
%   At each point, WENO4 considers a convex combination of two
%   quadratic interpolants, achieving a (bi,tri)-cubic interpolant in
%   smooth regions.  This uses the same stencil as a tri-cubic
%   interpolant.
%
%   w = weno4_interp(cpgrid, f, x)
%   Interpolates the data "f" on the grid given by "cpgrid" onto the
%   points "x".  There are certain assumptions about "x" and the grid,
%   namely that the band of the grid contains the stencil needed for
%   WENO4.
%
%   In the closest point method, the call would typically be:
%      w = weno4_interp(cpgrid, f, [cpx cpy cpz])
%
%   The dimension is determined from the number of columns of x.
%
%   "cpgrid" must contain fields cpgrid.x1d, .y1d, .band (and .z1d
%   in 3D).
%
%   For multiple calls with the same "x" (e.g., in the closest point
%   method) some data can be cached to improve performance:
%      wenoCache = weno4_interp(cpgrid, f, x, 'cache')
%      u1 = weno4_interp(wenoCache, f1)
%      u2 = weno4_interp(wenoCache, f2)
%      ...
%
%   The scheme implemented here is derived in [Macdonald & Ruuth
%   2008, Level Set Equations on Surfaces...].
%
%   TODO: dual-band support.

  if (nargin < 5)
    opt = [];
    opt.wenoeps = 1e-6;
  end

  if (nargin == 2)
    % have previous cached computations, stored in "cp"
    dim = cp.dim;
    if (dim == 2)
      w = weno4_interp2d_cached(cp, f, opt);
    else
      w = weno4_interp3d_cached(cp, f, opt);
    end
    return
  end

  if (nargin >= 4)
    if (makeCache)
      % noop
    elseif (strcmpi(makeCache, 'cache'))
      makeCache = true;
    else
      makeCache = false;
    end
  else
    makeCache = false;
  end

  [n1,dim] = size(x);

  if makeCache
    if dim == 2
      Cache = weno4_interp2d_makecache(cp, f, x);
    elseif dim == 3
      Cache = weno4_interp3d_makecache(cp, f, x);
    else
      error('dim not implemented');
    end
    w = Cache;
    return
  end

  % just do a single interpolation
  if (nargin == 3)
    if dim == 2
      w = weno4_interp2d(cp, f, x, opt);
    else
      w = weno4_interp3d(cp, f, x, opt);
    end
  end
end




function w = weno4_interp3d_cached(Cache, f, opt)
%WENO4_INTERP3D  nonlinear WENO interpolation 3D

  C = Cache.C;
  xi = Cache.xi;
  yi = Cache.yi;
  zi = Cache.zi;
  x = Cache.x;
  y = Cache.y;
  z = Cache.z;
  dx = Cache.dx;
  dy = Cache.dy;
  dz = Cache.dz;

  %tic
  u = {};
  v = {};
  for k=1:4
    for j=1:4
      u{j} = helper1d(C{1,j,k}*f, C{2,j,k}*f, C{3,j,k}*f, C{4,j,k}*f, xi, dx, x, opt);
    end
    v{k} = helper1d(u{1}, u{2}, u{3}, u{4}, yi, dy, y, opt);
  end
  w = helper1d(v{1}, v{2}, v{3}, v{4}, zi, dz, z, opt);
  %toc

end




function w = weno4_interp2d_cached(Cache, f, opt)
%WENO4_INTER23D  nonlinear WENO interpolation 3D

  C = Cache.C;
  xi = Cache.xi;
  yi = Cache.yi;
  x = Cache.x;
  y = Cache.y;
  dx = Cache.dx;
  dy = Cache.dy;

  %tic
  u = {};
  for j=1:4
    u{j} = helper1d(C{1,j}*f, C{2,j}*f, C{3,j}*f, C{4,j}*f, xi, dx, x, opt);
  end
  w = helper1d(u{1}, u{2}, u{3}, u{4}, yi, dy, y, opt);
  %toc
end




function u = helper1d(fim1, fi, fip1, fip2, xi, dx, x, opt)
%HELPER1D  A 1D WENO4 implementation

  WENOEPS = opt.wenoeps;  % the WENO parameter to prevent div-by-zero

  IS1 = ( 26*fip1.*fim1  -  52*fi.*fim1  -  76*fip1.*fi ...
          + 25*fip1.^2  +  64*fi.^2  +  13*fim1.^2 ) / 12;

  IS2 = ( 26*fip2.*fi  -  52*fip2.*fip1  -  76*fip1.*fi ...
          + 25*fi.^2  +  64*fip1.^2  +  13*fip2.^2 ) / 12;

  C1 = ((xi + 2*dx) - x) / 3*dx;
  C2 = (x - (xi - dx)) / 3*dx;
  alpha1 = C1 ./ (WENOEPS + IS1).^2;
  alpha2 = C2 ./ (WENOEPS + IS2).^2;
  w1 = alpha1 ./ (alpha1 + alpha2);
  w2 = alpha2 ./ (alpha1 + alpha2);

  p1 = fi  +  (x-xi).*(fip1 - fim1)/(2*dx)  +  (x-xi).^2 .* (fip1 - 2*fi + fim1)/(2*dx^2);
  p2 = fi  +  (x-xi).*(-fip2 + 4*fip1 - 3*fi)/(2*dx)  +  (x-xi).^2 .* (fip2 - 2*fip1 + fi)/(2*dx^2);

  u = w1.*p1 + w2.*p2;
end




function Cache = weno4_interp3d_makecache(cp, f, xyz)
%WENO4_INTERP3D_MAKE_CACHE  pre-computation for WENO interpolation 3D

  x1d = cp.x1d;
  y1d = cp.y1d;
  z1d = cp.z1d;
  Nx = length(x1d);
  Ny = length(y1d);
  Nz = length(z1d);
  dx = x1d(2) - x1d(1);
  dy = y1d(2) - y1d(1);
  dz = z1d(2) - z1d(1);

  relpt = [cp.x1d(1)  cp.y1d(1)  cp.z1d(1)];

  % determine the basepoint
  [ijk,X] = findGridInterpBasePt_vec(xyz, 3, relpt, [dx dy dz]);
  xi = X(:,1) + dx;
  yi = X(:,2) + dy;
  zi = X(:,3) + dz;
  ijk = ijk + 1;

  I = sub2ind([Ny Nx Nz], ijk(:,2), ijk(:,1), ijk(:,3));

  B = findInBand(I, cp.band, Nx*Ny*Nz);

  %tic
  [E W N S U D] = neighbourMatrices(cp, cp.band, cp.band);
  %toc

  %tic
  C = {};
  I = speye(size(W));
  T1 = {W,I,E,E*E};
  T2 = {S,I,N,N*N};
  T3 = {D,I,U,U*U};
  for i=1:4
    for j=1:4
      for k=1:4
        C{i,j,k} = B*T1{i}*T2{j}*T3{k};
      end
    end
  end
  %toc

  Cache.C = C;
  Cache.xi = xi;
  Cache.yi = yi;
  Cache.zi = zi;
  Cache.x = xyz(:,1);
  Cache.y = xyz(:,2);
  Cache.z = xyz(:,3);
  Cache.dx = dx;
  Cache.dy = dy;
  Cache.dz = dz;
  Cache.dim = 3;
end




function Cache = weno4_interp2d_makecache(cp, f, xyz)
%WENO4_INTERP2D_MAKE_CACHE  pre-computation for WENO interpolation 2D

  x1d = cp.x1d;
  y1d = cp.y1d;
  Nx = length(x1d);
  Ny = length(y1d);
  dx = x1d(2) - x1d(1);
  dy = y1d(2) - y1d(1);

  relpt = [cp.x1d(1)  cp.y1d(1)];

  % determine the basepoint
  [ijk,X] = findGridInterpBasePt_vec(xyz, 3, relpt, [dx dy]);
  xi = X(:,1) + dx;
  yi = X(:,2) + dy;
  ijk = ijk + 1;

  I = sub2ind([Ny Nx], ijk(:,2), ijk(:,1));

  B = findInBand(I, cp.band, Nx*Ny);

  [E W N S] = neighbourMatrices(cp, cp.band, cp.band);
  %toc

  %tic
  C = {};
  I = speye(size(W));
  T1 = {W,I,E,E*E};
  T2 = {S,I,N,N*N};
  %T3 = {D,I,U,U*U};
  for i=1:4
    for j=1:4
      C{i,j} = B*T1{i}*T2{j};
    end
  end
  %toc

  Cache.C = C;
  Cache.xi = xi;
  Cache.yi = yi;
  Cache.x = xyz(:,1);
  Cache.y = xyz(:,2);
  Cache.dx = dx;
  Cache.dy = dy;
  Cache.dim = 2;
end


function w = weno4_interp2d(cp, f, xyz, opt)
%WENO4_INTERP2D  nonlinear WENO interpolation 2D

  %tic
  x1d = cp.x1d;
  y1d = cp.y1d;
  Nx = length(x1d);
  Ny = length(y1d);
  dx = x1d(2) - x1d(1);
  dy = y1d(2) - y1d(1);

  relpt = [cp.x1d(1)  cp.y1d(1)];

  x = xyz(:,1);
  y = xyz(:,2);

  % determine the basepoint, roughly speaking this is "floor(xy)"
  % in terms of the grid
  [ijk,X] = findGridInterpBasePt_vec(xyz, 3, relpt, [dx dy]);
  % +1 here because the basepoint is actually the lowerleft corner
  % of the stencil and we want the "floor".
  xi = X(:,1) + dx;
  yi = X(:,2) + dy;
  ijk = ijk + 1;

  I = sub2ind([Ny Nx], ijk(:,2), ijk(:,1));
  B = findInBand(I, cp.band, Nx*Ny);

  [E W N S] = neighbourMatrices(cp, cp.band, cp.band);
  %preptime = toc;

  % some duplicated work because many interpolation points will have the
  % same basepoint
  %tic
  g = S*f;     u1 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = f;       u2 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = N*f;     u3 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = N*(N*f); u4 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);

  w = helper1d(u1, u2, u3, u4, yi, dy, y, opt);
  %wenotime = toc;
  %fprintf('weno4: preptime=%g, wenotime=%g\n', preptime, wenotime)
end



function w = weno4_interp3d(cp, f, xyz, opt)
%WENO4_INTERP3D  nonlinear WENO interpolation 3D

  %tic
  x1d = cp.x1d;
  y1d = cp.y1d;
  z1d = cp.z1d;
  Nx = length(x1d);
  Ny = length(y1d);
  Nz = length(z1d);
  dx = x1d(2) - x1d(1);
  dy = y1d(2) - y1d(1);
  dz = z1d(2) - z1d(1);

  relpt = [cp.x1d(1)  cp.y1d(1)  cp.z1d(1)];

  x = xyz(:,1);
  y = xyz(:,2);
  z = xyz(:,3);

  % determine the basepoint, roughly speaking this is "floor(xy)"
  % in terms of the grid
  [ijk,X] = findGridInterpBasePt_vec(xyz, 3, relpt, [dx dy dz]);
  % +1 here because the basepoint is actually the lowerleft corner
  % of the stencil and we want the "floor".
  xi = X(:,1) + dx;
  yi = X(:,2) + dy;
  zi = X(:,3) + dz;
  ijk = ijk + 1;

  I = sub2ind([Ny Nx Nz], ijk(:,2), ijk(:,1), ijk(:,3));
  B = findInBand(I, cp.band, Nx*Ny*Nz);

  [E W N S U D] = neighbourMatrices(cp, cp.band, cp.band);
  %preptime = toc;

  % some duplicated work because many interpolation points will have the
  % same basepoint

  %tic
  g = D*(S*f);     u1 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = D*f;         u2 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = D*(N*f);     u3 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = D*(N*(N*f)); u4 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  w1 = helper1d(u1, u2, u3, u4, yi, dy, y, opt);

  g = S*f;         u1 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = f;           u2 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = N*f;         u3 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = N*(N*f);     u4 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  w2 = helper1d(u1, u2, u3, u4, yi, dy, y, opt);

  g = U*(S*f);     u1 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = U*f;         u2 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = U*(N*f);     u3 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = U*(N*(N*f)); u4 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  w3 = helper1d(u1, u2, u3, u4, yi, dy, y, opt);

  g = U*(U*(S*f));     u1 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = U*(U*f);         u2 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = U*(U*(N*f));     u3 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  g = U*(U*(N*(N*f))); u4 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x, opt);
  w4 = helper1d(u1, u2, u3, u4, yi, dy, y, opt);

  w = helper1d(w1, w2, w3, w4, zi, dz, z, opt);
  %wenotime = toc;
  %fprintf('weno4: preptime=%g, wenotime=%g\n', preptime, wenotime)
end
