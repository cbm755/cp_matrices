function varargout = interp_matrix(g, p, other)
%INTERP_MATRIX  Return an interpolation matrix
%   E = interp_matrix(g, p)
%       Build a matrix which interpolates grid data on a grid defined
%       by the cpgrid object 'g' onto its closest points.
%       Interpolation is done using degree 'p' barycentric Lagrange
%       interpolation (defaults to 3 if omitted).  E will be a
%       'length(g.x)' by 'length(g.band)' sparse matrix
%
%   E = interp_matrix(g, {XI YI}, p)
%   E = interp_matrix(g, {XI YI ZI}, p)
%   E = interp_matrix(g, {XI YI ZI ... WI}, p)
%       As above but you can specify a cell array of column vectors
%       for the interpolation points.  Often used for the "Eplot"
%       matrix
%
%   [Ei, Ej, Es] = interp_matrix(...)
%       Here the entries of the matrix are returned as three vectors
%       (like in FEM).  This avoids some overhead of constructing the
%       matrix (and can be useful in determining computational bands).
%
%   If you need access to the full (non-banded) E matrix, use the
%   dimension-dependent code (e.g., inter2_matrix()).
%
%   If your grid object has an inverse bandmap (g.invbandmap),
%   it will be passed to the low level function.  TODO: this needs
%   documented somewhere.

  if ~isstruct(g)
    error('expected a cpgrid object');
  end

  if nargin == 1
    p = [];
  end

  if nargin == 3
    xi = p;
    p = other;
    if ~iscell(xi)
      error('expected cell array of interpolation points');
    end
  else
    xi = {};
  end

  if (isempty(p))
    p = 3;  % default interp degree
  end

  if iscell(g.x1d)
    if isfield(g, 'invbandmap')
      ibm = g.invbandmap;
    else
      ibm = [];
    end
    if isempty(xi)
      xi = g.cpx;
    end
    [varargout{1:nargout}] = ...
        interpn_matrix(g.x1d, xi, p, g.band, ibm);
  elseif g.dim == 2
    if isempty(xi)
      [varargout{1:nargout}] = ...
          interp2_matrix(g.x1d, g.y1d, g.cpx, g.cpy, p, g.band, false);
    else
      [varargout{1:nargout}] = ...
          interp2_matrix(g.x1d, g.y1d, xi{1}, xi{2}, p, g.band, false);
    end
  elseif g.dim == 3
    if isempty(xi)
      [varargout{1:nargout}] = ...
          interp3_matrix(g.x1d,g.y1d,g.z1d, g.cpx,g.cpy,g.cpz, p, ...
                         g.band, false);
    else
      [varargout{1:nargout}] = ...
          interp3_matrix(g.x1d,g.y1d,g.z1d, xi{1},xi{2},xi{3}, p, ...
                         g.band, false);
    end
  else
    error('this grid object cannot be handled');
  end

