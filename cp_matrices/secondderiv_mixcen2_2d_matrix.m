function [Dxyc] = secondderiv_mixcen2_2d_matrix(x,y, band1, band2, use_ndgrid)
%SECONDDERIV_MIXCEN2_2D_MATRIX  Build discrete second mixed partial
%   Matrices for 2nd-order centred differences, mixed partials, in 2D
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

  weights = [1  -1   1  -1] / (4*dx*dy);
  PTS = [ 1   1; ...
         -1   1; ...
         -1  -1; ...
          1  -1];
  Dxyc = helper_diff_matrix2d(x, y, band1, band2, weights, PTS, use_ndgrid);
