function L = laplacian_nd_matrix(xs, order, band1, band2)
%LAPLACIAN_ND_MATRIX  Build discrete Laplacian in n-D
% Supports 2nd- and 4th-order accurate
%
%   L = laplacian_nd_matrix(X, 2, band)
%   L = laplacian_nd_matrix(X, 2, band1, band2)
%   L = laplacian_nd_matrix(X, 4, band1, band2)
%
% Note only does ndgrid ordering.

  if (nargin < 3)
    band2 = band1;
  end
  if (nargin < 4)
    order = 2;
  end

  dim = length(xs);

  % just call the second deriv's
  D2c = secondderiv_cen2_nd_matrices(xs, band1, band2, order);

  % and add them up
  L = D2c{1};
  for n=2:dim
    L = L + D2c{n};
  end

