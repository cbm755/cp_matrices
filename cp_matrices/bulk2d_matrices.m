function [Dxx,Dyy,Dxc,Dyc,Dxb,Dyb,Dxf,Dyf,Dxyc] = bulk2d_matrices(x, y, use_ndgrid)
%BULK2D_MATRICES  Build discrete derivative matrices
%
%   To use ndgrid ordering pass "true" as the final argument

  if (nargin <= 2)
    use_ndgrid = false;
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);

  nx = length(x);
  ny = length(y);

  %% build 1D operators
  [Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(nx, dx);
  [Iy,D1yy,D1yc,D1yb,D1yf] = diff_matrices1d(ny, dy);

  % Use kronecker products to build 2D operators
  if (use_ndgrid)
    Dxx = kron(Iy, D1xx);
    Dyy = kron(D1yy, Ix);
    % laplacian
    %L = Dxx + Dyy;

    Dxyc = kron(D1yc, D1xc);

    Dxc = kron(Iy, D1xc);
    Dyc = kron(D1yc, Ix);

    Dxb = kron(Iy, D1xb);
    Dyb = kron(D1yb, Ix);

    Dxf = kron(Iy, D1xf);
    Dyf = kron(D1yf, Ix);
  else % meshgrid ordering
    Dxx = kron(D1xx, Iy);
    Dyy = kron(Ix, D1yy);

    Dxyc = kron(D1xc, D1yc);

    Dxc = kron(D1xc, Iy);
    Dyc = kron(Ix, D1yc);

    Dxb = kron(D1xb, Iy);
    Dyb = kron(Ix, D1yb);

    Dxf = kron(D1xf, Iy);
    Dyf = kron(Ix, D1yf);
  end
