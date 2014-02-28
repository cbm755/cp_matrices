function [varargout] = secondderiv_matrices(g, also)
%SECONDDERIV_MATRICES  Build discrete second derivatives
%
%   [Dxx, Dyy] = seconddeiv_matrices(g)
%   [Dxx, Dyy, Dxy] = seconddeiv_matrices(g)
%      Builds second-order accurate centered differences.
%
%   [Dxx, Dyy, Dzz] = secondderiv_matrices(g);
%   [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = secondderiv_matrices(g)
%      Works in 3D (here g uses a 3D meshgrid).
%
%   In general dimensions:
%     D = secondderiv_matrices(g);
%     [Dxx,Dyy,...,Dzz] = secondderiv_matrices(g);
%   n-dimenionsional version, here g contains an ndgrid() and D =
%   {Dxx,Dyy,...,Dzz} is a cell array of length n.
%   You can also get the mixed partials in n dimensions:
%     [Dxx,Dyy,...,Dzz] = secondderiv_matrices(g, 'nd_mixed');
%   TODO: implement this
%
%   TODO: band1, band2?
%

  if nargin < 2
    also = false;
  elseif strcmp(lower(also), 'nd_mixed')
    also = true;
  else
    also = false;
  end



  % detect general dimension (don't use g.dim b/c e.g., ndgrid in 3D)
  if iscell(g.x1d)  % arbitrary dimension, ndgird

    if isfield(g, 'invbandmap')
      ibm = g.invbandmap;
    else
      ibm = [];
    end

    if ~also  % just straight partials
      out = secondderiv_cen2_nd_matrices(g.x1d, g.band, g.band, ibm);
      if nargout == 1
        varargout{1} = out;
      else
        for i=1:g.dim
          varargout{i} = out{i};
        end
      end
    else % mixed partials
      error('TODO: not implemented');
    end

  elseif g.dim == 2  % 2D Meshgrid

    [varargout{1:2}] = ...
        secondderiv_cen2_2d_matrices(g.x1d,g.y1d,g.band,g.band,false);

    if nargout > 2
      varargout{3} = secondderiv_mixcen2_2d_matrix(g.x1d,g.y1d, ...
                                                   g.band,g.band,false);
    end

  elseif g.dim == 3  % 3D Meshgrid

    [varargout{1:3}] = secondderiv_cen2_3d_matrices( ...
        g.x1d, g.y1d, g.z1d, g.band, g.band, false);
    if nargout > 3
      [varargout{4:6}] = secondderiv_mixcen2_3d_matrices( ...
          g.x1d, g.y1d, g.z1d, g.band, g.band, false);
    end

  end

