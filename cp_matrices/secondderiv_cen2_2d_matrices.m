function [Dxxc, Dyyc] = secondderiv_cen2_2d_matrices(x,y, band1, band2, use_ndgrid)
%SECONDDERIV_CEN2_2D_MATRICES  Build discrete second derivatives
%   Matrices for 2nd-order centred differences in 2D
%
%   To use ndgrid ordering pass "true" as the final argument

  if (nargin <= 4)
    use_ndgrid = false;
  end
  if (nargin <= 3)
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

  weights = [1  -2  1] / dx^2;
  PTS = [-1   0; ...
          0   0; ...
          1   0];
  Dxxc = helper_diff_matrix2d(x, y, band1, band2, weights, PTS, use_ndgrid);

  weights = [1  -2  1] / dy^2;
  PTS = [ 0  -1; ...
          0   0; ...
          0   1];
  Dyyc = helper_diff_matrix2d(x, y, band1, band2, weights, PTS, use_ndgrid);

