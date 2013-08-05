function [Db, Df] = ...
  firstderiv_upw1_nd_matrices(xs, band1, band2, invbandmap)
%FIRSTDERIV_UPW1_ND_MATRICES  1st-order upwinded first derivatives
%   [Db, Df] = firstderiv_upw1_nd_matrices(X, band)
%   [Db, Df] = firstderiv_upw1_nd_matrices(X, band1, band2)
%
%   Db is a cell array of backward differences in each coordinate
%   direction.  Df is cell array of forward differences.
%
%   Note only does ndgrid ordering.

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

  Db = {};
  Df = {};

  for n=1:dim
    weights = [-1  1] / ddx(n);
    PTS = zeros(2, dim);
    PTS(1, n) = -1;
    Db{n} = helper_diff_matrixnd(NN, band1, band2, weights, PTS, invbandmap);

    weights = [-1  1] / ddx(n);
    PTS = zeros(2, dim);
    PTS(2, n) = 1;
    Df{n} = helper_diff_matrixnd(NN, band1, band2, weights, PTS, invbandmap);
  end

