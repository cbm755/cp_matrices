function [Dxx,Dyy,Dzz, Dxc,Dyc,Dzc, Dxb,Dyb,Dzb, Dxf,Dyf,Dzf] ...
      = bulk3d_matrices(x, y, z, use_ndgrid)
%BULK2D_MATRICES  Build discrete derivative matrices
%
%   To use ndgrid ordering pass "true" as the final argument

  if (nargin <= 3)
    use_ndgrid = false;
  end

  if use_ndgrid
    error('not implemented');
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);

  nx = length(x);
  ny = length(y);
  nz = length(z);

  %% build 1D operators
  [Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(nx, dx);
  [Iy,D1yy,D1yc,D1yb,D1yf] = diff_matrices1d(ny, dy);
  [Iz,D1zz,D1zc,D1zb,D1zf] = diff_matrices1d(nz, dz);



  % Use kronecker products to build 3D operators
  Dxx = kron(Iz, kron(D1xx, Iy));
  Dyy = kron(Iz, kron(Ix, D1yy));
  Dzz = kron(kron(D1zz, Ix), Iy);

  % laplacian
  %L = Dxx + Dyy + Dzz;

  Dxc = kron(Iz, kron(D1xc, Iy));
  Dyc = kron(Iz, kron(Ix, D1yc));
  Dzc = kron((kron(D1zc), Ix), Iy);

  warning('not tested/implemented yet');

  Dxb = 1;
  Dyb = 1;
  Dzb = 1;

  Dxf = 1;
  Dyf = 1;
  Dzf = 1;