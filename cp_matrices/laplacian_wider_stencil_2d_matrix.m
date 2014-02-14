function L = laplacian_wider_stencil_2d_matrix(x,y, order, alpha, beta, band1, band2, use_ndgrid, use_loop)
%LAPLACIAN_2D_MATRIX  Build a 2D discrete Laplacian
%   ORDER: 2 for 2nd-order
%   Does no error checking up the equispaced nature of x,y,z
%
%   TODO: explain the roles of band1 and band2
%
%   TODO: currently assumes dx=dy, but should be an easy fix
%
%   TODO: issue with two bands: currently need extra padding in
%   wherever we call this from: should fix this
%
%   To use ndgrid ordering pass "true" as the final argument

  if (nargin <= 8)
    use_loop = false;
  end
  if (nargin <= 7)
    use_ndgrid = false;
  end
  if (nargin <= 6)
    band2 = band1;
  end

  % input checking
  [temp1, temp2] = size(x);
  if ~(  (ndims(x) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('x must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(y);
  if ~(  (ndims(y) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('y must be a vector, not e.g., meshgrid output');
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  if ~assertAlmostEqual(dx, dy, 100*eps)
    error('this particular routine requires dx == dy');
  end
  %ddx = [dx  dy];
  %dim = length(ddx);
  %Nx = round( (x(end)-x(1)) / dx ) + 1;
  %Ny = round( (y(end)-y(1)) / dy ) + 1;
  %ptL = [x(1) y(1)];
  %ptH = [x(end) y(end)];


  if (order == 2)
    weights = [ -4*alpha/dx^2-2*beta/dx^2 ...
                alpha*[1 1 1 1] / dx^2 ...
                beta*[1 1 1 1] / (2*dx^2) ];
    PTS = [ 0   0; ...
            1   0; ...
           -1   0; ...
            0   1; ...
            0  -1; ...
            1   1; ...
           -1   1; ...
            1  -1; ...
           -1  -1];
  else
    error(['order ' num2str(order) ' not implemented']);
  end

  if (use_loop)
    L = helper_diff_matrix2d_oldloop(x, y, band1, band2, weights, PTS, use_ndgrid);
  else
    L = helper_diff_matrix2d(x, y, band1, band2, weights, PTS, use_ndgrid);
  end


