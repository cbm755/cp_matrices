function L = laplacian_wider_stencil_3d_matrix(x,y,z, order, alpha, beta, gamma, band1, band2, use_ndgrid, use_loop)
%LAPLACIAN_3D_MATRIX  Build a 3D discrete Laplacian
%   ORDER: 2 for 2nd -order
%   Does no error checking up the equispaced nature of x,y,z
%
%   L = laplacian_3d_matrix(x,y,z, order, band)
%      'L' is a discrete laplacian over a grid.
%      'order' can be 2 or 4.
%      x,y,z are 1d vectors which form a meshgrid, 'band' is a
%      subset of this meshgrid given as linear indices.
%      TODO: this will have a dirichlet boundary condition at the
%      edge of the band, at least for order 2.
%
%   L = laplacian_3d_matrix(x,y,z, order, band1, band2)
%      dual-banded version of above, a rectangular matrix.
%      TODO: explain this
%
%   To use ndgrid ordering pass "true" as an extra argument.
%
%   Pass "true" as a further argument to use the (slower)
%   looping-based code.

  if (nargin <= 10)
    use_loop = false;
  end
  if (nargin <= 9)
    use_ndgrid = false;
  end
  if (nargin <= 8)
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
  [temp1, temp2] = size(z);
  if ~(  (ndims(z) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('z must be a vector, not e.g., meshgrid output');
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);
  if ~assertAlmostEqual([dx dx], [dy dz], 100*eps)
    error('this particular routine requires dx == dy == dz');
  end
  %ddx = [dx  dy  dz];
  %dim = length(ddx);
  %Nx = round( (x(end)-x(1)) / dx ) + 1;
  %Ny = round( (y(end)-y(1)) / dy ) + 1;
  %Nz = round( (z(end)-z(1)) / dz ) + 1;
  %ptL = [x(1) y(1) z(1)];
  %ptH = [x(end) y(end) z(end)];

  if (order == 2)
    weights = [ -6 * ( alpha/dx^2 + beta/(3*dx^2) + gamma/(2*dx^2) ) ...
                alpha * ones(1,6) / dx^2 ...
                beta * 6/8 * ones(1,8) / (3*dx^2) ...
                gamma * 6/12 * ones(1,12) /(2*dx^2) ];
    PTS = [ 0   0   0; ...
% alpha: center of faces
            1   0   0; ...
           -1   0   0; ...
            0   1   0; ...
            0  -1   0; ...
            0   0   1; ...
            0   0  -1; ...
% beta:  corner points 
            1   1   1; ...
            1   1  -1; ...
            1  -1   1; ...
           -1   1   1; ...
           -1  -1   1; ...
           -1   1  -1; ...
            1  -1  -1; ...
           -1  -1  -1; ...
% gamma: center of edges
            0   1  -1; ...
            0   1   1; ...
            0  -1   1; ...
            0  -1  -1; ...
            1   0  -1; ...
            1   0   1; ...
           -1   0   1; ...
           -1   0  -1; ...
            1  -1   0; ...
            1   1   0; ...
           -1   1   0; ...
           -1  -1   0 ];

  else
    error(['order ' num2str(order) ' not implemented']);
  end

  if (use_loop)
    L = helper_diff_matrix3d_oldloop(x, y, z, band1, band2, weights, PTS, use_ndgrid);
  else
    L = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);
  end
