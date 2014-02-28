function [varargout] = firstderiv_matrices(g, which)
%FIRSTDERIV_MATRICES  Build discrete first derivatives
%
%   [Dxc,Dyc] = firstderiv_matrices(g)
%   [Dxc,Dyc] = firstderiv_matrices(g, 'centered')
%      Builds second-order accurate centered differences.
%
%   [Dxb,Dyb] = firstderiv_matrices(g, 'backward')
%      Builds 1st-order accurate backward difference matrices.
%
%   [Dxf,Dyf] = firstderiv_matrices(g, 'forward')
%      Builds 1st-order accurate forward difference matrices.
%
%   [Dx,Dy,Dz] = firstderiv_matrices(g, ...);
%      Works in 3D (here g uses a 3D meshgrid).
%
%   [D] = firstderiv_matrices(g, ...);
%   [Dx,Dy,...,Dw] = firstderiv_matrices(g, ...);
%      n-dimenionsional version, here g contains an ndgrid() and D =
%      {Dx,Dy,...,Dw} is a cell array of length n.
%
%
%   TODO: band1, band2?
%

  if nargin < 2
    which = 'centered';
  end

  % detect general dimension (don't use g.dim b/c e.g., ndgrid in 3D)
  if iscell(g.x1d)
    if isfield(g, 'invbandmap')
      ibm = g.invbandmap;
    else
      ibm = [];
    end
  end

  switch lower(which)
    case 'centered'
      if iscell(g.x1d)
        out = firstderiv_cen2_nd_matrices(g.x1d, g.band, g.band, ibm);
      elseif g.dim == 2
        [varargout{1:2}] = ...
            firstderiv_cen2_2d_matrices(g.x1d, g.y1d, g.band, g.band, false);
      elseif g.dim == 3
        [varargout{1:3}] = ...
            firstderiv_cen2_3d_matrices(g.x1d, g.y1d, g.z1d, g.band, g.band, false);
      end

    case 'forward'
      % TODO: For historical reasons, this implementation generates both
      % forward and backward matrices and simply disgards the other
      % one.  Inefficient, and should be easy to fix.
      if iscell(g.x1d)
        [Db, Df] = firstderiv_upw1_nd_matrices(g.x1d, g.band, g.band, ibm);
        %varargout{1} = Df;
        out = Df;
        % TODO: expanded output here?
      elseif g.dim == 2
        [Dxb,Dxf,Dyb,Dyf] = ...
            firstderiv_upw1_2d_matrices(g.x1d, g.y1d, g.band, g.band, false);
        varargout{1} = Dxf;
        varargout{2} = Dyf;
      elseif g.dim == 3
        [Dxb,Dxf,Dyb,Dyf,Dzb,Dzf] = ...
            firstderiv_upw1_3d_matrices(g.x1d,g.y1d,g.z1d,g.band,g.band,false);
        varargout{1} = Dxf;
        varargout{2} = Dyf;
        varargout{3} = Dzf;
      end

    case 'backward'
      if iscell(g.x1d)
        [Db, Df] = firstderiv_upw1_nd_matrices(g.x1d, g.band, g.band, ibm);
        out = Db;
        %varargout{1} = Db;
      elseif g.dim == 2
        [Dxb,Dxf,Dyb,Dyf] = ...
            firstderiv_upw1_2d_matrices(g.x1d, g.y1d, g.band, g.band, false);
        varargout{1} = Dxb;
        varargout{2} = Dyb;
      elseif g.dim == 3
        [Dxb,Dxf,Dyb,Dyf,Dzb,Dzf] = ...
            firstderiv_upw1_3d_matrices(g.x1d,g.y1d,g.z1d,g.band,g.band,false);
        varargout{1} = Dxb;
        varargout{2} = Dyb;
        varargout{3} = Dzb;
      end

    otherwise
      error('invalid grid input');
  end

  % Cell array or list output for nd case
  if iscell(g.x1d)
    if nargout == 1
      varargout{1} = out;
    else
      for i=1:g.dim
        varargout{i} = out{i};
      end
    end
  end