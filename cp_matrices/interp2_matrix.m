function [E,Ej,Es] = interp2_matrix(x, y, xi, yi, p, band, use_ndgrid)
%INTERP2_MATRIX  Return a 2D interpolation matrix
%   E = INTERP2_MATRIX(X, Y, XI, YI, P)
%   Build a matrix which interpolates grid data on a grid defined by
%   the product of the lists X and Y onto the points specified by
%   the lists XI and YI.  Interpolation is done using degree P
%   barycentric Lagrange interpolation.  E will be a length(XI) times
%   length(X)*length(Y) sparse matrix.  If P is omitted or set to
%   [] it defaults to 3.
%
%   E = INTERP2_MATRIX(X, Y, XI, YI, P, BAND)
%   BAND is a list of linear indices into a (possibly fictious) 2D
%   array of points constructed with meshgrid.  Here the columns of E
%   that are not in BAND are discarded.  This is done by first
%   constructing E as above.  E will be a length(XI) times
%   length(BAND) sparse matrix.  If BAND is [], no banding is done.
%
%   [Ei,Ej,Es] = INTERP2_MATRIX(X, Y, XI, YI, P)
%   [Ei,Ej,Es] = INTERP2_MATRIX(X, Y, XI, YI, P, BAND)
%   Here the entries of the matrix are returned as three vectors
%   (like in FEM).  This is efficient and avoids the overhead of
%   constructing the matrix.  If BAND is passed or not determines
%   the column space of the result (i.e., effects Ej).
%
%   INTERP2_MATRIX(X, Y, XI, YI, P, BAND, true)
%      Uses ndgrid instead of meshgrid ordering.
%
%   Notes:
%      This is faster replacement for INTERP2_MATRIX_OLDLOOP.
%      Inputs X and Y must be equispaced but this is not checked.

%   TODO:
%      interp2, interp3 and interpn have a lot of common code,
%      particularly the Ei,Ej,Es to matrix bits which are dimension
%      independent.  Should refactor.


  % input checking
  if ~isvector(x) || ~isvector(y)
    error('x and y must be vectors, not e.g., meshgrid output');
  end
  if ~(  (ndims(xi) == 2) && (size(xi,2) == 1)  )
    error('xi must be a column vector');
  end
  if ~(  (ndims(yi) == 2) && (size(yi,2) == 1)  )
    error('yi must be a column vector');
  end

  if (nargin == 4)
    p = [];
    makeBanded = false;
    use_ndgrid = false;
  elseif (nargin == 5)
    makeBanded = false;
    use_ndgrid = false;
  elseif (nargin == 6)
    if isempty(band) makeBanded = false; else makeBanded = true; end
    use_ndgrid = false;
  elseif (nargin == 7)
    if isempty(band) makeBanded = false; else makeBanded = true; end
  else
    error('unexpected inputs');
  end

  if (isempty(p))
    p = 3;
  end

  if (nargout > 1)
    makeListOutput = true;
  else
    makeListOutput = false;
  end

  T1 = cputime();
  dx = x(2)-x(1);   Nx = length(x);
  dy = y(2)-y(1);   Ny = length(y);
  ddx = [dx  dy];
  ptL = [x(1) y(1)];
  M = Nx*Ny;

  if (M > 1e15)
    error('too big to use doubles as indicies: implement int64 indexing')
  end

  dim = 2;
  N = p+1;
  EXTSTENSZ = N^dim;

  %tic
  Ei = repmat((1:length(xi))',1,EXTSTENSZ);
  Ej = zeros(size(Ei));
  weights = zeros(size(Ei));
  % todo: integers seem slower(!), although use less memory.
  %Ei = repmat(uint32((1:length(xi))'),1,EXTSTENSZ);
  %Ej = zeros(size(Ei), 'uint32');
  %weights = zeros(size(Ei), 'double');
  %toc

  %tic
  % this used to be a call to buildInterpWeights but now most of
  % that is done here
  [Ibpt, Xgrid] = findGridInterpBasePt_vec({xi yi}, p, ptL, ddx);
  xw = LagrangeWeights1D_vec(Xgrid{1}, xi, ddx(1), N);
  yw = LagrangeWeights1D_vec(Xgrid{2}, yi, ddx(2), N);
  %toc

  %tic
  % this is a good order for memory access: ijk just counts up
  for i=1:N
    for j=1:N
      gi = (Ibpt{1} + i - 1);
      gj = (Ibpt{2} + j - 1);
      ijk = sub2ind([N,N], j, i);
      weights(:,ijk) = xw(:,i) .* yw(:,j);

      if (use_ndgrid)
        Ej(:,ijk) = sub2ind([Nx,Ny], gi, gj);
      else
        %Ej(:,ijk) = sub2ind([Ny,Nx], gj, gi);
        Ej(:,ijk) = (gi-1)*Ny + gj;
      end
    end
  end
  %toc
  T1 = cputime() - T1;
  %fprintf('done new Ei,Ej,weights, total time: %g\n', T1);

  % TODO: is there any advantage to keeping Ei as matrices?  Then each
  % column corresponds to the same point in the stencil...
  if ~makeListOutput
    %tic
    E = sparse(Ei(:), Ej(:), weights(:), length(xi), M);
    % contruct the transpose instead: faster and less memory but
    % makeBanded is very slow below on the transpose
    %E = sparse(Ej(:), Ei(:), weights(:), M, length(xi));
    %T2 = toc;
    %fprintf('call to "sparse" time: %g\n', toc);

    if (makeBanded)
      nnzEfull = nnz(E);
      E = E(:,band);
      %E = E(band,:);  % very slow for the transpose

      if nnz(E) < nnzEfull
        % Sanity check: the columns outside of band should all be
        % zero.  TODO: should be an error?
        warning('non-zero coefficients discarded by banding');
      end
    end
  end

  if (makeListOutput)
    E = Ei(:);   % first output is called E
    Ej = Ej(:);
    Es = weights(:);

    if (makeBanded)
      invbandmap = make_invbandmap(M, band);
      Ej = invbandmap(Ej);
    end
  end
end
