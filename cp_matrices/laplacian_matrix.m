function L = laplacian_matrix(g, order)
%LAPLACIAN_MATRIX  Build discrete Laplacian on grid
%   L = laplacian_matrix(g)
%      Builds a Laplacian over the grid 'g'.
%   L = laplacian_matrix(g, order)
%      order can be '2' or '4' and defaults to 2 if omitted.
%
%   TODO: band1, band2?

  if nargin < 2
    order = 2;
  end

  if iscell(g.x1d)
    if isfield(g, 'invbandmap')
      ibm = g.invbandmap;
    else
      ibm = [];
    end
    L = laplacian_nd_matrix(g.x1d, order, g.band, g.band, ibm);
  elseif g.dim == 2
    L = laplacian_2d_matrix(g.x1d, g.y1d, order, g.band, g.band, false);
  elseif g.dim == 3
    L = laplacian_3d_matrix(g.x1d, g.y1d, g.z1d, order, g.band, g.band, false);
  end

