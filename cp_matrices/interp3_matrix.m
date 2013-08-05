function [E,Ej,Es] = interp3_matrix(x, y, z, xi, yi, zi, p, band, use_ndgrid)
%INTERP3_MATRIX  Return a 3D interpolation matrix
%   E = INTERP3_MATRIX(X,Y,Z, XI,YI,ZI, P)
%   Build a matrix which interpolates grid data on a grid defined by
%   the product of the lists X, Y and Z onto the points specified by
%   the lists XI, YI, ZI.  Interpolation is done using degree P
%   barycentric Lagrange interpolation.  E will be a length(XI) times
%   length(X)*length(Y)*length(Z) sparse matrix.  If P is omitted or
%   set to [] it defaults to 3.
%
%   E = INTERP3_MATRIX(X,Y,Z, XI,YI,ZI, P, BAND)
%   BAND is a list of linear indices into a (possibly fictious) 3D
%   array of points constructed with meshgrid.  Here the columns of E
%   that are not in BAND are discarded.  This is done by first
%   constructing E as above.  E will be a length(XI) times
%   length(BAND) sparse matrix.  If BAND is [], no banding is done.
%
%   [Ei,Ej,Es] = INTERP3_MATRIX(X,Y,Z, XI,YI,ZI, P)
%   [Ei,Ej,Es] = INTERP3_MATRIX(X,Y,Z, XI,YI,ZI, P, BAND)
%   Here the entries of the matrix are returned as three vectors
%   (like in FEM).  This is efficient and avoids the overhead of
%   constructing the matrix.  If BAND is passed or not determines
%   the column space of the result (i.e., effects Ej).
%
%   INTERP3_MATRIX(X, Y, XI, YI, P, BAND, true)
%      Uses ndgrid instead of meshgrid ordering.
%
%   Notes:
%      This is faster replacement for INTERP3_MATRIX_OLDLOOP.
%      Inputs X, Y, Z must be equispaced but this is not checked.


  % input checking
  if (~isvector(x)) || (~isvector(y)) || (~isvector(z))
    error('x, y and z must be vectors, not e.g., meshgrid output');
  end
  if ~(  (ndims(xi) == 2) && (size(xi,2) == 1)  )
    error('xi must be a column vector');
  end
  if ~(  (ndims(yi) == 2) && (size(yi,2) == 1)  )
    error('yi must be a column vector');
  end
  if ~(  (ndims(zi) == 2) && (size(zi,2) == 1)  )
    error('zi must be a column vector');
  end

  if (nargin == 6)
    p = [];
    makeBanded = false;
    use_ndgrid = false;
  elseif (nargin == 7)
    makeBanded = false;
    use_ndgrid = false;
  elseif (nargin == 8)
    if isempty(band) makeBanded = false; else makeBanded = true; end
    use_ndgrid = false;
  elseif (nargin == 9)
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
  dz = z(2)-z(1);   Nz = length(z);
  ddx = [dx  dy  dz];
  ptL = [x(1) y(1) z(1)];
  M = Nx*Ny*Nz;

  if (M > 1e15)
    error('too big to use doubles as indicies: implement int64 indexing')
  end

  dim = 3;
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
  [Ibpt, Xgrid] = findGridInterpBasePt_vec({xi yi zi}, p, ptL, ddx);
  xw = LagrangeWeights1D_vec(Xgrid{1}, xi, ddx(1), N);
  yw = LagrangeWeights1D_vec(Xgrid{2}, yi, ddx(2), N);
  zw = LagrangeWeights1D_vec(Xgrid{3}, zi, ddx(3), N);
  %toc

  %tic
  % this is a good order for memory access: ijk just counts up
  for k=1:N
    for i=1:N
      for j=1:N
        gi = (Ibpt{1} + i - 1);
        gj = (Ibpt{2} + j - 1);
        gk = (Ibpt{3} + k - 1);
        ijk = sub2ind([N,N,N], j, i, k);
        weights(:,ijk) = xw(:,i) .* yw(:,j) .* zw(:,k);

        if (use_ndgrid)
          Ej(:,ijk) = sub2ind([Nx,Ny,Nz], gi, gj, gk);
          %Ej(:,ijk) = (gk-1)*(Nx*Ny) + (gj-1)*Nx + gi;
        else
          % all these do the same, but last one is fastest.  Although sub2ind
          % presumably has safety checks...
          %ind = (gk-1)*(Nx*Ny) + (gi-1)*Ny + gj;
          %ind = sub2ind([Ny,Nx,Nz], gj, gi, gk);
          %ind = round((gk-1)*(Nx*Ny) + (gi-1)*(Ny) + gj-1 + 1);
          Ej(:,ijk) = (gk-1)*(Nx*Ny) + (gi-1)*Ny + gj;
        end
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
    %T2 = toc;
    %fprintf('call to "sparse" time: %g\n', toc);

    if (makeBanded)
      nnzEfull = nnz(E);
      E = E(:,band);

      if nnz(E) < nnzEfull
        % sanity check: the columns outside of band should all be
        % zero.  TODO: should be an error?
        warning('non-zero coefficients discarded by banding');
      end
    end
  end

  if (1==0)
    disp('[testing] get back components:');
    tic; [I,J,V] = find(Es); toc

    disp('call "sparse" on smaller system:');
    tic; Es2 = sparse(I, J, V, length(xi), length(band)); toc
    E-Es2
  end

  if (1==0)
    % TODO: do this outside in another function
    disp('banding 2: this way finds innerband');
    %tic
    %innerband = unique(Ej(:));
    %toc
    tic
    [innerband,I,J] = unique(Ej(:));
    toc
    tic
    Es3 = sparse(Ei(:), J, weights(:), length(xi),length(innerband));
    toc
    %keyboard
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
