function [L,Lj,Ls] = laplacian_3d_matrix_tempcomp(x,y,z, order, band1, band2, use_ndgrid, use_loop)
%LAPLACIAN_3D_MATRIX  Build a 3D discrete Laplacian
%   ORDER: 2 or 4 for 2nd or 4th-order
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

  if (nargin <= 7)
    use_loop = false;
  end
  if (nargin <= 6)
    use_ndgrid = false;
  end
  if (nargin <= 5)
    band2 = band1;
  end

  if (nargout > 1)
    makeListOutput = true;
  else
    makeListOutput = false;
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
    weights = [-6 1 1 1 1 1 1] / dx^2;
    PTS = [ 0   0   0; ...
            1   0   0; ...
           -1   0   0; ...
            0   1   0; ...
            0  -1   0; ...
            0   0   1; ...
            0   0  -1];
  elseif (order == 4)
    weights = [-15.0/2.0 ...
               (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
               (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
               (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
              ] / dx^2;
    PTS = [  0   0   0; ...
            -2   0   0; ...
            -1   0   0; ...
             1   0   0; ...
             2   0   0; ...
             0  -2   0; ...
             0  -1   0; ...
             0   1   0; ...
             0   2   0; ...
             0   0  -2; ...
             0   0  -1; ...
             0   0   1; ...
             0   0   2];
  else
    error(['order ' num2str(order) ' not implemented']);
  end

  warning('experimental code');
  if (use_loop)
    % TODO: time to drop this?  Is it in unit tests though?
    L = helper_diff_matrix3d_oldloop(x, y, z, band1, band2, weights, PTS, use_ndgrid);
  else
    if makeListOutput
      warning('should fix helper instead');
      [L,Lj,Ls] = helper_diff_matrix3d_tempcomp(...
          x, y, z, band1, band2, weights, PTS, use_ndgrid);
    else
      L = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);
    end
  end
