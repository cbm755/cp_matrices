function [Dxb,Dxf, Dyb,Dyf, Dzb,Dzf] = ...
  firstderiv_upw1_3d_matrices(x,y,z, band1, band2, use_ndgrid)
%FIRSTDERIV_UPW1_3D_MATRICES  1st-order upwinded first derivatives
%   [Dxb,Dxf, Dyb,Dyf, Dzb,Dzf] = ...
%      firstderiv_upw1_3d_matrices(x,y,z, band)
%   Dxb is backward differences in x direction,
%   Dxf is forward differences in the x direction, etc.
%
%   [...] = firstderiv_upw1_3d_matrices(x,y,z, band1, band2)
%   Dual banded.
%
%   [...] = firstderiv_upw1_3d_matrices(x,y,z, band1, band2, true)
%   Pass 'true' as final argument to use ndgrid ordering.

  if (nargin <= 5)
    use_ndgrid = false;
  end
  if (nargin <= 4)
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
  %dim = length(ddx);

  weights = [-1  1] / dx;
  PTS = [-1   0   0; ...
          0   0   0];
  Dxb = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [-1  1] / dx;
  PTS = [ 0   0   0; ...
          1   0   0];
  Dxf = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [-1  1] / dy;
  PTS = [ 0  -1   0; ...
          0   0   0];
  Dyb = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [-1  1] / dy;
  PTS = [ 0   0   0; ...
          0   1   0];
  Dyf = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [-1  1] / dz;
  PTS = [ 0   0  -1; ...
          0   0   0];
  Dzb = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [-1  1] / dz;
  PTS = [ 0   0   0; ...
          0   0   1];
  Dzf = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);
