function [Dxx,Dyy,Dzz, Dxc,Dyc,Dzc, Dxb,Dyb,Dzb, Dxf,Dyf,Dzf, ...
          Dxyc,Dxzc,Dyzc] = bulk3d_matrices(x, y, z, use_ndgrid, BC)
%BULK3D_MATRICES  Build discrete derivative matrices
%
%   Usage:
%   [Dxx,Dyy,Dzz, Dxc,Dyc,Dzc, Dxb,Dyb,Dzb, Dxf,Dyf,Dzf, Dxyc,Dxzc,Dyzc] = ...
%      bulk3d_matrices(x, y, z, use_ndgrid, BC)
%
%   By default, this assumes meshgrid ordering and periodic BCs.
%
%   To use ndgrid ordering pass "true" as the fourth argument.
%
%   Select different BCs with the last argument (see "help
%   diff_matrices1d").

  if (nargin <= 3)
    use_ndgrid = false;
  end
  if (nargin <= 4)
    BC = 'p';
  end
  if isempty(use_ndgrid)
    use_ndgrid = false;
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);

  nx = length(x);
  ny = length(y);
  nz = length(z);

  %% build 1D operators
  [Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(nx, dx, BC);
  [Iy,D1yy,D1yc,D1yb,D1yf] = diff_matrices1d(ny, dy, BC);
  [Iz,D1zz,D1zc,D1zb,D1zf] = diff_matrices1d(nz, dz, BC);



  % Use kronecker products to build 3D operators
  if (use_ndgrid)
    % note those with kron(Iy,Ix) are indep of meshgrid/ndgrid ordering
    Dxx = kron(Iz, kron(Iy, D1xx));
    Dyy = kron(Iz, kron(D1yy, Ix));
    Dzz = kron(D1zz, kron(Iy, Ix));%

    Dxc = kron(Iz, kron(Iy, D1xc));
    Dyc = kron(Iz, kron(D1yc, Ix));
    Dzc = kron(D1zc, kron(Iy, Ix));%

    Dxb = kron(Iz, kron(Iy, D1xb));
    Dyb = kron(Iz, kron(D1yb, Ix));
    Dzb = kron(D1zb, kron(Iy, Ix));%

    Dxf = kron(Iz, kron(Iy, D1xf));
    Dyf = kron(Iz, kron(D1yf, Ix));
    Dzf = kron(D1zf, kron(Iy, Ix));%

    Dxyc = kron(Iz, kron(D1yc, D1xc));
    Dxzc = kron(D1zc, kron(Iy, D1xc));
    Dyzc = kron(D1zc, kron(D1yc, Ix));

  else
    Dxx = kron(Iz, kron(D1xx, Iy));
    Dyy = kron(Iz, kron(Ix, D1yy));
    Dzz = kron(D1zz, kron(Ix, Iy));
    % laplacian
    %L = Dxx + Dyy + Dzz;

    Dxc = kron(Iz, kron(D1xc, Iy));
    Dyc = kron(Iz, kron(Ix, D1yc));
    Dzc = kron(D1zc, kron(Ix, Iy));

    Dxb = kron(Iz, kron(D1xb, Iy));
    Dyb = kron(Iz, kron(Ix, D1yb));
    Dzb = kron(D1zb, kron(Ix, Iy));

    Dxf = kron(Iz, kron(D1xf, Iy));
    Dyf = kron(Iz, kron(Ix, D1yf));
    Dzf = kron(D1zf, kron(Ix, Iy));

    Dxyc = kron(Iz, kron(D1xc, D1yc));
    Dxzc = kron(D1zc, kron(D1xc, Iy));
    Dyzc = kron(D1zc, kron(Ix, D1yc));
  end