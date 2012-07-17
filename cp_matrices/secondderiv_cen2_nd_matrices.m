function [D2c] = secondderiv_cen2_nd_matrices(xs, band1, band2)
%SECONDDERIV_CEN2_ND_MATRICES  Build discrete first derivatives
% Matrices for 2nd-derivatives which are 2nd-order centered
% differences in n dimensions.
%
%   D2c = secondderiv_cen2_nd_matrices(X, band)
%   D2c = secondderiv_cen2_nd_matrices(X, band1, band2)
%
% D2c is a cell array of centered differences in each coordinate
% direction.
%
% TODO: support a call like:
%   [Dxc,Dyc,...Dwc] = firstderiv_cen2_nd_matrices(X, band1, band2)
%
% Note only does ndgrid ordering.

  if (nargin < 3)
    band2 = band1;
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
    weights = [1  -2   1] / ddx(n)^2;
    PTS = zeros(3, dim);
    PTS(1, n) = -1;
    PTS(3, n) = 1;
    D2c{n} = helper_diff_matrixnd(NN, band1, band2, weights, PTS);
  end

