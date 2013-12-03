function L = laplacian_3d_matrix(x,y,z, order, band1, band2, use_ndgrid, use_loop)
%LAPLACIAN_3D_MATRIX  Build a 3D discrete Laplacian
%
%   DEPRECATED?: you may want laplacian_matrix()
%
%   L = laplacian_3d_matrix(x,y,z, order, band)
%      'L' is a discrete laplacian over a grid.
%      'order' can be 2 or 4.
%      x,y,z are 1D vectors which form a meshgrid, 'band' is a
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
%   Pass "true" as a further argument to use the (older, slower)
%   looping-based code.
%
%   Currently if dx is within 10*macheps of dy and dz, this assumes
%   dx==dy==dz, otherwise it uses a more general formula.  The logic
%   is to have marginally less rounding error in the common dx==dy==dz
%   case: is it worth it?
%
%   Does no error checking up the equispaced nature of x,y,z.

  if (nargin <= 7)
    use_loop = false;
  end
  if (nargin <= 6)
    use_ndgrid = false;
  end
  if (nargin <= 5)
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
  if assertAlmostEqual([dx dx], [dy dz], 10*eps)
    dxequal = 1;
  else
    dxequal = 0;
  end

  if (order == 2)
    if dxequal
      weights = [-6 1 1 1 1 1 1] / dx^2;
    else
      weights = [ -2/dx^2 - 2/dy^2 - 2/dz^2  ...
                  [1 1]/dx^2   [1 1]/dy^2   [1 1]/dz^2 ];
    end
    PTS = [ 0   0   0; ...
            1   0   0; ...
           -1   0   0; ...
            0   1   0; ...
            0  -1   0; ...
            0   0   1; ...
            0   0  -1];
  elseif (order == 4)
    if dxequal
      weights = [-15.0/2.0 ...
                 (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
                 (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
                 (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
                ] / dx^2;
    else
      weights = [ -5/(2*dx^2) - 5/(2*dy^2) - 5/(2*dz^2)  ...
                  [-1/12  4/3  4/3  -1/12] / dx^2  ...
                  [-1/12  4/3  4/3  -1/12] / dy^2  ...
                  [-1/12  4/3  4/3  -1/12] / dz^2 ];
    end
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

  if (use_loop)
    L = helper_diff_matrix3d_oldloop(x, y, z, band1, band2, weights, PTS, use_ndgrid);
  else
    L = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);
  end

