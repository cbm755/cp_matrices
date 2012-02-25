function [E,Ej,Es] = interp2_matrix(x, y, xi, yi, p, band)
%INTERP2_MATRIX  Return a 2D interpolation matrix
%   E = INTERP2_MATRIX(X, Y, XI, YI, P)
%   Build a matrix which interpolates grid data on a grid defined by
%   the product of the lists X and Y onto the points specified by
%   the lists XI and YI.  Interpolation is done using degree P
%   barycentric Lagrange interpolation.  E will be a length(XI) times
%   length(X)*length(Y) sparse matrix.
%
%   E = INTERP2_MATRIX(X, Y, XI, YI, P, BAND)
%   BAND is a list of linear indices into a (possibly fictious) 2D
%   array of points constructed with meshgrid.  Here the columns of E
%   that are not in BAND are discarded.  This is done by first
%   constructing E as above.  E will be a length(XI) times
%   length(BAND) sparse matrix.
%
%   [Ei,Ej,Es] = INTERP2_MATRIX(X, Y, XI, YI, P)
%   [Ei,Ej,Es] = INTERP2_MATRIX(X, Y, XI, YI, P, BAND)
%   Here the entries of the matrix are returned as three vectors
%   (like in FEM).  This is efficient and avoids the overhead of
%   constructing the matrix.  If BAND is passed or not determines
%   the column space of the result (i.e., effects Ej).
%   (TODO: with BAND currently not implemented).
%
%   Does no error checking up the equispaced nature of x and y
%
%   Notes: this is faster replacement for INTERP2_MATRIX_OLDLOOP.

  % input checking
  [temp1, temp2] = size(x);
  if ~(  (ndims(x) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('x must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(y);
  if ~(  (ndims(y) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('y must be a vector, not e.g., meshgrid output');
  end
  if ~(  (ndims(xi) == 2) && (size(xi,2) == 1)  )
    error('xi must be a column vector');
  end
  if ~(  (ndims(yi) == 2) && (size(yi,2) == 1)  )
    error('yi must be a column vector');
  end

  if (nargin == 4)
    p = 3
    makeBanded = false;
  elseif (nargin == 5)
    makeBanded = false;
  elseif (nargin == 6)
    makeBanded = true;
  else
    error('unexpected inputs');
  end

  if (nargout > 1)
    makeListOutput = true;
  else
    makeListOutput = false;
  end

  if makeBanded && makeListOutput
    error('currently cannot make both Banded and Ei,Ej,Es output');
  end

  T = tic;
  dx = x(2)-x(1);   Nx = length(x);
  dy = y(2)-y(1);   Ny = length(y);
  ddx = [dx  dy];
  ptL = [x(1) y(1)];

  if (Nx * Ny > 1e15)
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
  [Ibpt, Xgrid] = findGridInterpBasePt_vec([xi yi], p, ptL, ddx);
  xw = LagrangeWeights1D_vec(Xgrid(:,1), xi, ddx(1), N);
  yw = LagrangeWeights1D_vec(Xgrid(:,2), yi, ddx(2), N);
  %toc

  %tic
  % this is a good order for memory access: ijk just counts up
  for i=1:N
    for j=1:N
      gi = (Ibpt(:,1) + i - 1);
      gj = (Ibpt(:,2) + j - 1);
      ijk = sub2ind([N,N], j, i);
      weights(:,ijk) = xw(:,i) .* yw(:,j);

      %Ej(:,ijk) = sub2ind([Ny,Nx], gj, gi);
      Ej(:,ijk) = (gi-1)*Ny + gj;
    end
  end
  %toc
  T1 = toc(T);
  %fprintf('done new Ei,Ej,weights, total time: %g\n', T1);


  % TODO: is there any advantage to keeping Ei as matrices?  Then each
  % column corresponds to the same point in the stencil...
  if ~makeListOutput
    tic
    E = sparse(Ei(:), Ej(:), weights(:), length(xi), Nx*Ny);
    T2 = toc;
    %fprintf('call to "sparse" time: %g\n', toc);
  end

  if (makeBanded)
    %disp('band the large matrix:');
    if (1==1)
      %tic
      E = E(:,band);
      %toc
    else
      % sanity check: the columns outside of band should all be zero
      tic
      Esparse = E(:,band);
      Eout = E(:,setdiff(1:(Nx*Ny),band));
      if (nnz(Eout) > 0)
        nnz(Eout)
        warning('Lost some non-zero coefficients (from outside the innerband)');
      end
      E = Esparse;
      toc
    end
  end

  if (makeListOutput)
    E = Ei(:);   % first output is called E
    Ej = Ej(:);
    Es = weights(:);
  end
end
