function w = weno4_interp_caching(cp, f, x, makeCache)
%WENO4_INTERP  nonlinear WENO interpolation in 2D/3D
%   At each point, WENO4 considers a convex combination of two
%   quadratic interpolants, achieving a (bi,tri)-cubic interpolant in
%   smooth regions.  This uses the same stencil as a tri-cubic
%   interpolant.
%
%   u = weno4_interp(cpgrid, f, x)
%   Interpolates the data "f" on the grid given by "cpgrid" onto the
%   points "x".  There are certain assumptions about "x" and the grid,
%   namely that the band of the grid contains the stencil needed for
%   WENO4.  In the closest point method, x might be [cpx cpy cpz].
%
%   The dimension is determined from the number of columns of x.
%
%   "cpgrid" must contain fields cpgrid.x1d, .y1d, .band (and .z1d
%   in 3D)
%
%   For a large system in 3D, this will be very slow.  But assuming
%   interpolation onto the same points will be repeated (as in the
%   closest point method), significant savings are possible by
%   cachine:
%   wenoCache = weno4_interp(cpgrid, f, x, 'cache')
%   u1 = weno4_interp(wenoCache, f1)
%   u2 = weno4_interp(wenoCache, f2)
%   ...
%
%   The scheme implemented here is derived in [Macdonald & Ruuth
%   2008, Level Set Equations on Surfaces...].
%
%   TODO: support nonvector relpt
%   TODO: support calling without a "cpgrid"?
%   TODO: dual-band support.

  if (nargin == 2)
    % have previous cached computations, stored in "cp"
    dim = cp.dim;
    if (dim == 2)
      w = weno4_interp2d(cp, f);
    else
      w = weno4_interp3d(cp, f);
    end
    return
  end

  if (nargin == 4)
    disp('weno4: bulding cache');
  end

  [n1,dim] = size(x);
  if dim == 2
    Cache = weno4_interp2d_makecache(cp, f, x);
  elseif dim == 3
    Cache = weno4_interp3d_makecache(cp, f, x);
  else
    error('dim not implemented');
  end

  if (nargin == 3)
    if dim == 2
      w = weno4_interp2d(f, Cache);
    else
      w = weno4_interp3d(f, Cache);
    end
  else
    w = Cache;
  end
end




function w = weno4_interp3d(Cache, f)
%WENO4_INTERP3D  nonlinear WENO interpolation 3D

  C = Cache.C;
  xi = Cache.xi;
  yi = Cache.yi;
  zi = Cache.zi;
  x = Cache.x;
  y = Cache.y;
  z = Cache.z;
  dx = Cache.dx;

  tic
  u = {};
  v = {};
  for k=1:4
    for j=1:4
      u{j} = helper1d(C{1,j,k}*f, C{2,j,k}*f, C{3,j,k}*f, C{4,j,k}*f, xi, dx, x);
    end
    v{k} = helper1d(u{1}, u{2}, u{3}, u{4}, yi, dx, y);
  end
  w = helper1d(v{1}, v{2}, v{3}, v{4}, zi, dx, z);
  toc

end




function w = weno4_interp2d(Cache, f)
%WENO4_INTER23D  nonlinear WENO interpolation 3D

  tic
  C = Cache.C;
  xi = Cache.xi;
  yi = Cache.yi;
  x = Cache.x;
  y = Cache.y;
  dx = Cache.dx;
  toc

  tic
  u = {};
  for j=1:4
    u{j} = helper1d(C{1,j}*f, C{2,j}*f, C{3,j}*f, C{4,j}*f, xi, dx, x);
  end
  w = helper1d(u{1}, u{2}, u{3}, u{4}, yi, dx, y);
  toc
end




function u = helper1d(fim1, fi, fip1, fip2, xi, dx, x)
%HELPER1D  A 1D WENO4 implementation

  WENOEPS = 1e-6;  % the WENO parameter to prevent div-by-zero

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

  dx = x1d(2) - x1d(1);  % assumed constant and same in x,y,z

  relpt = cp.x1d(1);  % TODO

  x = xyz(:,1);
  y = xyz(:,2);
  z = xyz(:,3);

  tic
  % determine the basepoint
  [ijk,X] = findGridInterpBasePt(xyz, 3, relpt, dx);
  xi = X(:,1) + dx;
  yi = X(:,2) + dx;
  zi = X(:,3) + dx;
  ijk = ijk + 1;

  I = sub2ind([Ny Nx Nz], ijk(:,2), ijk(:,1), ijk(:,3));
  B = findInBand(I, cp.band, Nx*Ny*Nz);
  toc

  tic
  [E W N S U D] = neighbourMatrices(cp, cp.band, cp.band);
  toc

  tic
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
  toc

  Cache.C = C;
  Cache.xi = xi;
  Cache.yi = yi;
  Cache.zi = zi;
  Cache.x = x;
  Cache.y = y;
  Cache.z = z;
  Cache.dx = dx;
  Cache.dim = 3;
end




function Cache = weno4_interp2d_makecache(cp, f, xy)
%WENO4_INTERP2D_MAKE_CACHE  pre-computation for WENO interpolation 2D

  x1d = cp.x1d;
  y1d = cp.y1d;
  Nx = length(x1d);
  Ny = length(y1d);

  dx = x1d(2) - x1d(1);  % assumed constant and same in x,y

  relpt = cp.x1d(1);  % TODO

  x = xy(:,1);
  y = xy(:,2);

  tic
  % determine the basepoint
  [ijk,X] = findGridInterpBasePt(xy, 3, relpt, dx);
  xi = X(:,1) + dx;
  yi = X(:,2) + dx;
  ijk = ijk + 1;

  I = sub2ind([Ny Nx], ijk(:,2), ijk(:,1));

  B = findInBand(I, cp.band, Nx*Ny*Nz);

  [E W N S] = neighbourMatrices(cp, cp.band, cp.band);
  toc

  tic
  C = {};
  I = speye(size(W));
  T1 = {W,I,E,E*E};
  T2 = {S,I,N,N*N};
  %T3 = {D,I,U,U*U};
  for i=1:4
    for j=1:4
      C{i,j,k} = B*T1{i}*T2{j};
    end
  end
  toc

  Cache.C = C;
  Cache.xi = xi;
  Cache.yi = yi;
  Cache.x = x;
  Cache.y = y;
  Cache.dx = dx;
  Cache.dim = 2;
end

