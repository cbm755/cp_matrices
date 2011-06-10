function [Dxyc, Dxzc, Dyzc] = firstderiv_mixcen2_3d_matrices(x,y,z, band1, band2, use_ndgrid)
%FIRSTDERIV_MIXCEN2_MATRICES3D  Build discrete first derivatives
%   Matrices for 2nd-order centred differences in 3D (mixed partials)
%
%   To use ndgrid ordering pass "true" as the final argument

  if (nargin <= 5)
    use_ndgrid = false;
  end
  if (nargin <= 4)
    band2 = band1;
  end


  % input checking
  [temp1, temp2] = size(x);
  if ~(  (ndims(x) == 2) & (temp1 == 1 | temp2 == 1)  )
    error('x must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(y);
  if ~(  (ndims(y) == 2) & (temp1 == 1 | temp2 == 1)  )
    error('y must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(z);
  if ~(  (ndims(z) == 2) & (temp1 == 1 | temp2 == 1)  )
    error('z must be a vector, not e.g., meshgrid output');
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);

  weights = [1  -1   1  -1] / (4*dx*dy);
  PTS = [ 1   1   0; ...
         -1   1   0; ...
         -1  -1   0; ...
          1  -1   0];
  Dxyc = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [1  -1   1  -1] / (4*dx*dz);
  PTS = [ 1   0   1; ...
         -1   0   1; ...
         -1   0  -1; ...
          1   0  -1];
  Dxzc = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);

  weights = [1  -1   1  -1] / (4*dx*dy);
  PTS = [0   1   1; ...
         0  -1   1; ...
         0  -1  -1; ...
         0   1  -1];
  Dyzc = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, use_ndgrid);
