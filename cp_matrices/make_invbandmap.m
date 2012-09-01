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
%      invbandmap = make_invbandmap({x1d, y1d, ..., w1d}, band);
%      invbandmap = make_invbandmap(NN, band);
%      invbandmap = make_invbandmap(cpgrid);
%
%   x1d, y1d, etc, are the 1-D vectors (which specify the underlying
%   grid via a meshgrid or ndgrid---a tensor product).  NN is the
%   vector [length(x1d) length(y1d) ... length(w1d)]
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
    NN = [length(x1d) length(y1d)];

  elseif (nargin == 4)
    dim = 3;
    x1d = varargin{1};
    y1d = varargin{2};
    z1d = varargin{3};
    band = varargin{4};
    NN = [length(x1d) length(y1d) length(z1d)];

  elseif (nargin == 2)  % n-D w/ cell array or a list NN
    if (iscell(varargin{1}))
      dim = 'n';  % just not 2 or 3
      X1d = varargin{1};
      dim = length(X1d);
      NN = [];
      for d=1:dim
        NN(d) = length(X1d{d});
      end
    else
      % list of lengths provided directly
      NN = varargin{1};
    end
    band = varargin{2};

  elseif (nargin == 1)  % cpgrid object input
    cpgrid = varargin{1};
    band = cpgrid.band;
    if (iscell(cpgrid.x1d))
      dim = length(cpgrid.x1d);
      NN = [];
      for d=1:dim
        NN(d) = length(cpgrid.x1d{d});
      end
    else
      dim = cpgrid.dim;
      if (dim == 2)
        NN = [length(cpgrid.x1d) length(cpgrid.y1d)];
      elseif (dim == 3)
        NN = [length(cpgrid.x1d) length(cpgrid.y1d) length(cpgrid.z1d)];
      else
        error('must use n-D form');
      end
    end
  end

  invbandmap = make_invbandmap_private(NN, band);
end



function invbandmap = make_invbandmap_private(NN, band)
%MAKE_INVBANDMAP_PRIVATE  Given a vector of the lengths of x1d, etc
%and the band, build in the inverse band map

  M = 1;
  for i = 1:length(NN)
    M = M * NN(i);
  end

  % this is (I, J, a_{IJ}).  J=1 expands to ones(1,length(band))
  sparselin2bandmap = sparse(band, 1, 1:length(band), M,1);
  % the full() is here because the result is "big sparse"
  invbandmap = @(i) full(sparselin2bandmap(i));
end
