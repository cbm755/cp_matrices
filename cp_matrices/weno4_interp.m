function w = weno4_interp(cp, f, x)
%WENO4_INTERP  nonlinear WENO interpolation in 2D/3D
%   At each point, WENO4 considers a convex combination of two
%   quadratic interpolants, achieving a (bi,tri)-cubic interpolant in
%   smooth regions.  This uses the same stencil as a tri-cubic
%   interpolant.
%
%   w = weno4_interp(cpgrid, f, x)
%   Interpolates the data "f" on the grid given by "cpgrid" onto
%   the points "x".  There are certain assumptions about "x" and
%   the grid, namely that the band of the grid contains the stencil
%   needed for WENO4.
%
%   In the closest point method, the call would typically be:
%   w = weno4_interp(cpgrid, f, [cpx cpy cpz])
%
%   The dimension is determined from the number of columns of x.
%
%   "cpgrid" must contain fields cpgrid.x1d, .y1d, .band (and .z1d
%   in 3D)
%
%   For multiple calls wih the same "x" (e.g., in the closest point
%   method) WENO4_INTERP_CACHING should be about twice as fast.  On
%   the other hand, if you just do one interpolation, than this is
%   probably faster.
%
%   The scheme implemented here is derived in [Macdonald & Ruuth
%   2008, Level Set Equations on Surfaces...].
%
%   TODO: support vector relpt
%   TODO: support calling without a "cpgrid"?
%   TODO: dual-band support.
%   TODO: assumes dx = dy = dx

  [n1,dim] = size(x);
  if dim == 2
    w = weno4_interp2d(cp, f, x);
  elseif dim == 3
    w = weno4_interp3d(cp, f, x);
  else
    error('dim not implemented');
  end
end


function w = weno4_interp2d(cp, f, xy)
%WENO4_INTERP2D  nonlinear WENO interpolation 2D

  tic
  x1d = cp.x1d;
  y1d = cp.y1d;
  Nx = length(x1d);
  Ny = length(y1d);

  dx = x1d(2) - x1d(1);  % assumed constant and same in x,y

  relpt = cp.x1d(1);  % TODO: support nonvector relpt

  x = xy(:,1);
  y = xy(:,2);

  % determine the basepoint, roughly speaking this is "floor(xy)"
  % in terms of the grid
  [ij,X] = findGridInterpBasePt(xy, 3, relpt, dx);
  % +1 here because the basepoint is actually the lowerleft corner
  % of the stencil and we want the "floor".
  xi = X(:,1) + dx;
  yi = X(:,2) + dx;
  ij = ij + 1;

  I = sub2ind([Ny Nx], ij(:,2), ij(:,1));
  B = findInBand(I, cp.band, Nx*Ny);

  [E W N S] = neighbourMatrices(cp, cp.band, cp.band);
  preptime = toc;

  % some duplicated work because many interpolation points will have the
  % same basepoint
  tic
  g = S*f;     u1 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = f;       u2 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = N*f;     u3 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = N*(N*f); u4 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);

  w = helper1d(u1, u2, u3, u4, yi, dx, y);
  wenotime = toc;
  fprintf('weno4: preptime=%g, wenotime=%g\n', preptime, wenotime)
end



function w = weno4_interp3d(cp, f, xyz)
%WENO4_INTERP3D  nonlinear WENO interpolation 3D

  tic
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

  % determine the basepoint, roughly speaking this is "floor(xy)"
  % in terms of the grid
  [ijk,X] = findGridInterpBasePt(xyz, 3, relpt, dx);
  % +1 here because the basepoint is actually the lowerleft corner
  % of the stencil and we want the "floor".
  xi = X(:,1) + dx;
  yi = X(:,2) + dx;
  zi = X(:,3) + dx;
  ijk = ijk + 1;

  I = sub2ind([Ny Nx Nz], ijk(:,2), ijk(:,1), ijk(:,3));
  B = findInBand(I, cp.band, Nx*Ny*Nz);

  [E W N S U D] = neighbourMatrices(cp, cp.band, cp.band);
  preptime = toc;

  % some duplicated work because many interpolation points will have the
  % same basepoint

  % could cache some of these matrix products
  tic
  g = D*(S*f);     u1 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = D*f;         u2 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = D*(N*f);     u3 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = D*(N*(N*f)); u4 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  w1 = helper1d(u1, u2, u3, u4, yi, dx, y);

  g = S*f;         u1 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = f;           u2 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = N*f;         u3 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = N*(N*f);     u4 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  w2 = helper1d(u1, u2, u3, u4, yi, dx, y);

  g = U*(S*f);     u1 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = U*f;         u2 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = U*(N*f);     u3 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = U*(N*(N*f)); u4 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  w3 = helper1d(u1, u2, u3, u4, yi, dx, y);

  g = U*(U*(S*f));     u1 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = U*(U*f);         u2 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = U*(U*(N*f));     u3 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  g = U*(U*(N*(N*f))); u4 = helper1d(B*(W*g), B*g, B*(E*g), B*(E*(E*g)), xi, dx, x);
  w4 = helper1d(u1, u2, u3, u4, yi, dx, y);

  w = helper1d(w1, w2, w3, w4, zi, dx, z);
  wenotime = toc;

  fprintf('weno4: preptime=%g, wenotime=%g\n', preptime, wenotime)
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