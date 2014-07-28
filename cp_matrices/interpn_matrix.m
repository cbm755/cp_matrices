function [Ei, Ej, Es] = interpn_matrix(xs, xi, p, band, invbandmap)
%INTERPN_MATRIX  Return a n-D interpolation matrix
%
%   You probably want the "interp_matrix" frontend.
%
%   E = INTERPN_MATRIX(X, XI, P)
%   E = INTERPN_MATRIX({x y z ... w}, {xi yi zi ... wi}, P)
%   Build a matrix which interpolates grid data on a grid defined by
%   the product of the lists X onto the points specified by lists XI.
%   Interpolation is done using degree P barycentric Lagrange
%   interpolation.  E will be a 'length(xi)' by M sparse matrix
%   where M is the product of the lengths of X, Y, Z, ..., W.
%
%   E = INTERPN_MATRIX({x y z ... w}, [xi yi zi ... wi], P)
%   Alternatively, second input can be a matrix with columns
%   specifying the interpolation points (somewhat deprecated).
%
%   E = INTERPN_MATRIX(X, XI, P, BAND)
%   BAND is a list of linear indices into a (possibly fictious) n-D
%   array of points constructed with *ndgrid*.  Here the columns of E
%   that are not in BAND are discarded.  E will be a 'size(XI,1)' by
%   'length(BAND)' sparse matrix.
%
%   [Ei,Ej,Es] = INTERPN_MATRIX(...)
%   Here the entries of the matrix are returned as three vectors
%   (like in FEM).  This is efficient and avoids the overhead of
%   constructing the matrix.  If BAND is passed or not determines
%   the column space of the result (i.e., effects Ej).
%
%   E = INTERPN_MATRIX(X, XI, P, BAND, INVBANDMAP)
%   If you have an inverse band map, you can pass it as the last
%   argument to avoid reconstructing it
%
%   This code assumes the grid is equispaced but it does not do
%   error checking on this.

  if ~iscell(xs)
    error('expected a cell array of {x1d,y1d,...}');
  end

  if (nargin < 5)
    invbandmap = [];
  end
  if (nargin < 4)
    band = [];
  end
  if (nargin < 3)
    p = [];
  end

  if (nargout > 1)
    makeListOutput = true;
  else
    makeListOutput = false;
  end

  if (isempty(p))
    p = 3;  % default interp degree
  end

  if isempty(band)
    makeBanded = false;
  else
    makeBanded = true;
  end


  T1 = cputime();
  dim = length(xs);
  Ns = zeros(1, dim);
  ddx = zeros(1, dim);
  ptL = zeros(1, dim);
  for d=1:dim
    Ns(d) = length(xs{d});
    ddx(d) = xs{d}(2)-xs{d}(1);
    ptL(d) = xs{d}(1);
  end
  M = prod(Ns);
  if (M > 1e15)
    error('too big to use doubles as indicies: implement int64 indexing')
  end

  if ~iscell(xi)
    warning('do you want a cell array for interp pts {xi,yi,...}?');
    xi2 = cell(1,d);
    for d=1:dim
      xi2{d} = xi(:,d);
    end
    xi = xi2;
  end

  if makeBanded && isempty(invbandmap)
    invbandmap = make_invbandmap(M, band);
  end

  Nsten = p+1;
  EXTSTENSZ = (Nsten)^dim;

  Ni = length(xi{1}(:));

  Ei = repmat((1:Ni)',1,EXTSTENSZ);
  Ej = zeros(size(Ei));
  Es = zeros(size(Ei));

  [Ibpt, Xgrid] = findGridInterpBasePt_vec(xi, p, ptL, ddx);
  xw = {};
  for d=1:dim
    xw{d} = LagrangeWeights1D_vec(Xgrid{d}, xi{d}, ddx(d), Nsten);
  end

  NN = Nsten*ones(1,dim);
  for s=1:prod(NN)
    [ii{1:dim}] = ind2sub(NN, s);
    %weights(:,s) = xw(:,i) .* yw(:,j) .* zw(:,k);
    temp = xw{1}(:,ii{1});
    for d=2:dim
      temp = temp .* xw{d}(:,ii{d});
    end
    Es(:,s) = temp;

    for d=1:dim
      gi{d} = (Ibpt{d} + ii{d} - 1);
    end

    Ej(:,s) = sub2ind(Ns, gi{:});
  end
  T1 = cputime() - T1;
  %fprintf('done new Ei,Ej,weights, total time: %g\n', T1);


  % TODO: old memory hungry implementation: this is faster but not
  % practical in high dimension.  Deprecated?  (see helper_diff_nd)
  if (1==0)
    if ~makeListOutput
      E = sparse(Ei(:), Ej(:), Es(:), Ni, M);
      if (makeBanded)
        nnzEfull = nnz(E);
        E = E(:,band);
        if nnz(E) < nnzEfull
          warning('non-zero coefficients discarded by banding');
        end
      end
      return
    end
  end


  if (1==0)
    % TODO: do this outside in another function
    disp('banding 2: this way finds innerband');
    [innerband,I,J] = unique(Ej(:));
    Es3 = sparse(Ei(:), J, weights(:), length(xi),length(innerband));
  end

  Ei = Ei(:);
  Ej = Ej(:);
  Es = Es(:);

  if (makeBanded)
    Ej = invbandmap(Ej);
    numcols = length(band);
  else
    numcols = M;
  end

  % make the matrix unless user wants list output
  if ~makeListOutput
    Ei = sparse(Ei, Ej, Es, Ni, numcols);
  end
end

