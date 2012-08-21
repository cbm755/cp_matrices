function invbandmap = make_invbandmap(varargin)
%MAKE_INVBANDMAP  Make a mapping from natural indices to the band
%   The "band" in the closest point method gives a mapping from any
%   point in the computational domain to the "natural" linear index in
%   the (theoretical) underlying grid.  Sometimes it is useful to be
%   able to invert this map.  This routine returns a function (via a
%   function handle) which can be used to perform the mapping of
%   linear indices to band points.
%
%      invbandmap = make_invbandmap(x1d, y1d, band);
%      invbandmap = make_invbandmap(x1d, y1d, z1d, band);
%      invbandmap = make_invbandmap({x1d, y1d, ..., z1d}, band);
%      invbandmap = make_invbandmap(cpgrid);
%
%   x1d, y1d, etc, are the 1-D vectors (which specify the underlying
%   grid via a meshgrid or ndgrid---a tensor product).
%
%   Examples:
%     >> invbandmap = make_invbandmap(cpgrid);
%     >> invbandmap(cpgrid.band(101))
%     ans = 101
%
%     >> invbandmap(1)
%     ans = 0
%     (assuming 1 is not a point in the band)
%
%   Notes:
%     A hash-table is a good choice to store this map.  But
%     Matlab's containers.Map does not seem to support vector
%     lookups (that is, querying an array of keys without a loop).


  if (nargin == 3)
    dim = 2;
    x1d = varargin{1};
    y1d = varargin{2};
    band = varargin{3};
  elseif (nargin == 4)
    dim = 3;
    x1d = varargin{1};
    y1d = varargin{2};
    z1d = varargin{3};
    band = varargin{4};
  elseif (nargin == 2)  % n-D
    dim = 'n';  % just not 2 or 3
    X1d = varargin{1};
    assert(iscell(X1d));
    band = varargin{2};
  elseif (nargin == 1)  % cpgrid object input
    cpgrid = varargin{1};
    band = cpgrid.band;
    if (iscell(cpgrid.x1d))
      dim = 'n';
      X1d = cpgrid.x1d;
    else
      dim = cpgrid.dim;
      if (dim == 2)
        x1d = cpgrid.x1d;
        y1d = cpgrid.y1d;
      elseif (dim == 3)
        x1d = cpgrid.x1d;
        y1d = cpgrid.y1d;
        z1d = cpgrid.z1d;
      else
        error('must use n-D form');
      end
    end
  end

  if (dim == 2)
    M = length(x1d) * length(y1d);
  elseif (dim == 3)
    M = length(x1d) * length(y1d) * length(z1d);
  else
    M = 1;
    for i = 1:length(X1d)
      M = M * length(X1d{i});
    end
  end

  % this is (I, J, a_{IJ}).  J=1 expands to ones(1,length(band))
  sparselin2bandmap = sparse(band, 1, 1:length(band), M,1);
  invbandmap = @(i) full(sparselin2bandmap(i));

