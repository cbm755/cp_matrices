function [Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x,y,z, band1, band2, use_ndgrid)
%FIRSTDERIV_CEN2_3D_MATRICES  2nd-order centered first derivatives
%   [Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x, y, z, band)
%   generates the matrices that perform central differencing in the x-
%   and y-directions respectively. Grid spacing is included in the
%   matrix.
%
%   [...] = firstderiv_cen2_3d_matrices(x, y, z, band1, band2)
%   performs the same, with an optional outer band.
%
%   [...] = firstderiv_cen2_3d_matrices(x,y,z,band1,band2,use_ndgrid)
%   with 'use_ndgrid' set as 'true' uses ndgrid ordering instead of
%   meshgrid ordering.


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

  weights = [-1./2  0  1./2] / dx;
  PTS = [-1   0   0; ...
          0   0   0; ...
          1   0   0];
  Dxc = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [-1./2  0  1./2] / dy;
  PTS = [ 0  -1   0; ...
          0   0   0; ...
          0   1   0];
  Dyc = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [-1./2  0  1./2] / dz;
  PTS = [ 0   0  -1; ...
          0   0   0; ...
          0   0   1];
  Dzc = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);
