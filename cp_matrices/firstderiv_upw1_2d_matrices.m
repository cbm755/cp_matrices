function [Dxb,Dxf, Dyb,Dyf] = firstderiv_upw1_2d_matrices(x,y, band1, band2, use_ndgrid)
%FIRSTDERIV_UPW1_2D_MATRICES  1st-order upwinded first derivatives
%   [Dxb,Dxf, Dyb,Dyf] = firstderiv_upw1_2d_matrices(x,y, band)
%   Dxb is backward differences in x direction,
%   Dxf is forward differences in the x direction, etc.
%
%   [...] = firstderiv_upw1_3d_matrices(x,y, band1, band2)
%   Dual banded.
%
%   [...] = firstderiv_upw1_3d_matrices(x,y, band1, band2, true)
%   Pass 'true' as final argument to use ndgrid ordering.

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

  weights = [-1  1] / dx;
  PTS = [-1   0; ...
          0   0];
  Dxb = helper_diff_matrix2d(x, y, band1, band2, weights, PTS, use_ndgrid);

  weights = [-1  1] / dx;
  PTS = [ 0   0; ...
          1   0];
  Dxf = helper_diff_matrix2d(x, y, band1, band2, weights, PTS, use_ndgrid);

  weights = [-1  1] / dy;
  PTS = [ 0  -1; ...
          0   0];
  Dyb = helper_diff_matrix2d(x, y, band1, band2, weights, PTS, use_ndgrid);

  weights = [-1  1] / dy;
  PTS = [ 0   0; ...
          0   1];
  Dyf = helper_diff_matrix2d(x, y, band1, band2, weights, PTS, use_ndgrid);
