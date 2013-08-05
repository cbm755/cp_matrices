function [Dc] = firstderiv_cen2_nd_matrices(xs, band1, band2, invbandmap)
%FIRSTDERIV_CEN2_ND_MATRICES  2nd-order centered first derivatives
%   Matrices for 1st-derivatives which are 2nd-order centered
%   differences in n dimensions.
%
%   Dc = firstderiv_cen2_nd_matrices(X, band)
%   Dc = firstderiv_cen2_nd_matrices(X, band1, band2)
%
%   Dc is a cell array of centered differences in each coordinate
%   direction.
%
%   Note only does ndgrid ordering.
%
%   TODO: support a call like:
%   [Dxc,Dyc,...Dwc] = firstderiv_cen2_nd_matrices(X, band1, band2)

  if (nargin < 3)
    band2 = band1;
  end
  if (nargin < 4)
    invbandmap = [];
  end

  % TODO: input checking

  dim = length(xs);
  Ns = zeros(1, dim);
  ddx = zeros(1, dim);
  for n=1:dim
    NN(n) = length(xs{n});
    ddx(n) = xs{n}(2)-xs{n}(1);
  end

  Dc = {};

  for n=1:dim
    weights = [-1./2  0  1./2] / ddx(n);
    PTS = zeros(3, dim);
    PTS(1, n) = -1;
    PTS(3, n) = 1;
    Dc{n} = helper_diff_matrixnd(NN, band1, band2, weights, PTS, invbandmap);
  end

