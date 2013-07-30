function L = laplacian_nd_matrix(xs, order, band1, band2, invbandmap)
%LAPLACIAN_ND_MATRIX  Build discrete Laplacian in n-D
% Supports 2nd- and 4th-order accurate
%
%   L = laplacian_nd_matrix(X, 2, band)
%   L = laplacian_nd_matrix(X, 2, band1, band2)
%   L = laplacian_nd_matrix(X, 4, band1, band2)
%   L = laplacian_nd_matrix(X, 2, band1, band2, invbandmap)
%
% Note only does ndgrid ordering.

  if (nargin < 4)
    band2 = band1;
  end
  if (nargin < 2)
    order = 2;
  end

  dim = length(xs);
  NN = zeros(1, dim);
  ddx = zeros(1, dim);
  for n=1:dim
    NN(n) = length(xs{n});
    ddx(n) = xs{n}(2)-xs{n}(1);
  end

  if (nargin < 5)
    % slight speedup to calculate this once and pass it around
    invbandmap = make_invbandmap(prod(NN), band2);
  end

  D2c = secondderiv_cen2_nd_matrices(xs, band1, band2, order, invbandmap);

  dxAllEqual = all(ddx(2:end) == ddx(1));

  % Add up all the second derivatives
  if dxAllEqual
    % TODO: ugly hack to reduce roundoff.  Better: modify the helper
    % routine to use dx=1, then add, then divide by dx^2.
    L = (ddx(1)^2)*D2c{1};
    for n=2:dim
      L = L + (ddx(n)^2)*D2c{n};
    end
    L = L / ddx(1)^2;
  else
    L = D2c{1};
    for n=2:dim
      L = L + D2c{n};
    end
  end


