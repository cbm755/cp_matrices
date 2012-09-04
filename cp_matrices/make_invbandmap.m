function invbandmap = make_invbandmap(M, band)
%MAKE_INVBANDMAP  Make a mapping from natural indices to the band
%   The "band" in the closest point method gives a mapping from any
%   point in the computational domain to the "natural" linear index in
%   the (theoretical) underlying grid.  Sometimes it is useful to be
%   able to invert this map.  This routine returns a function (via a
%   function handle) which can be used to perform the mapping of
%   linear indices to band points.
%
%      invbandmap = make_invbandmap(cpgrid);
%
%   Alternatively:
%
%      M = length(x1d)*length(y1d);
%      M = prod([length(x1d) ... length(z1d)];
%      invbandmap = make_invbandmap(M, band);
%
%   where M is the total number of points in the meshgrid.  Here x1d,
%   y1d, etc, are the 1-D vectors (which specify the underlying grid
%   via a tensor product with meshgrid or ndgrid).
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
%     Here we have used a sparse matrix instead.  It must be a
%     column matrix in Matlab because Matlab allocates memory for
%     each *column* in a sparse matrix (and the trick is avoiding
%     allocated anything of size M, which might be very large).


  if (nargin == 1)  % cpgrid object input
    cpgrid = M;
    band = cpgrid.band;
    if (iscell(cpgrid.x1d))
      dim = length(cpgrid.x1d);
      M = 1;
      for d=1:dim
        M = M * length(cpgrid.x1d{d});
      end
    else
      dim = cpgrid.dim;
      if (dim == 2)
        M = length(cpgrid.x1d)*length(cpgrid.y1d);
      elseif (dim == 3)
        M = length(cpgrid.x1d)*length(cpgrid.y1d)*length(cpgrid.z1d);
      else
        error('must use n-D (cell-array) form');
      end
    end
  end

  if ~(isscalar(M))
    % list of lengths provided
    M = prod(varargin{1});
  end

  % this is (I, J, a_{IJ}).  J=1 expands to ones(1,length(band))
  sparselin2bandmap = sparse(band, 1, 1:length(band), M,1);
  % the full() is here because the result is "big sparse"
  invbandmap = @(i) full(sparselin2bandmap(i));



  % Build a mapping from logical grid to the band, various ways to do
  % this.  Here are the approaches I tried:
  %switch 1
  %  case 1  % lowest memory usage (sparse column vector)
  %    logical2bandmap = sparse(band, 1, 1:length(band), M,1);
  %  case 2
  %    logical2bandmap = sparse(band, 1, 1:length(band), M,1);
  %    % transposing it will use much more memory (4*Nx*Ny*Nz bytes)
  %    % but applying the operator is faster
  %    logical2bandmap = logical2bandmap';
  %  case 3  % construct the transpose directly, similar to 2 I think
  %    logical2bandmap = sparse(1, band, 1:length(band), 1,M);
  %  case 4  % use a dense map: fastest but most memory
  %    logical2bandmap = zeros(M,1);
  %    logical2bandmap(band) = 1:length(band);
  %end
  %Ej = logical2bandmap(Ej);


