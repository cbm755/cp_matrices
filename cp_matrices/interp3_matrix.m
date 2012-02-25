function [E,Ej,Es] = interp3_matrix(x, y, z, xi, yi, zi, p, band)
%INTERP3_MATRIX  Return a 3D interpolation matrix
%   E = INTERP3_MATRIX(X,Y,Z, XI,YI,ZI, P)
%   Build a matrix which interpolates grid data on a grid defined by
%   the product of the lists X, Y and Z onto the points specified by
%   the lists XI, YI, ZI.  Interpolation is done using degree P
%   barycentric Lagrange interpolation.  E will be a length(XI) times
%   length(X)*length(Y)*length(Z) sparse matrix.
%
%   E = INTERP3_MATRIX(X,Y,Z, XI,YI,ZI, P, BAND)
%   BAND is a list of linear indices into a (possibly fictious) 3D
%   array of points constructed with meshgrid.  Here the columns of E
%   that are not in BAND are discarded.  This is done by first
%   constructing E as above.  E will be a length(XI) times
%   length(BAND) sparse matrix.
%
%   [Ei,Ej,Es] = INTERP3_MATRIX(X,Y,Z, XI,YI,ZI, P)
%   [Ei,Ej,Es] = INTERP3_MATRIX(X,Y,Z, XI,YI,ZI, P, BAND)
%   Here the entries of the matrix are returned as three vectors
%   (like in FEM).  This is efficient and avoids the overhead of
%   constructing the matrix.  If BAND is passed or not determines
%   the column space of the result (i.e., effects Ej).
%   (TODO: with BAND currently not implemented).
%
%   Does no error checking up the equispaced nature of x,y,z
%
%   Notes: this is faster replacement for INTERP3_MATRIX and
%   INTERP3_MATRIX_BAND.


  % input checking
  [temp1, temp2] = size(x);
  if ~(  (ndims(x) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('x must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(y);
  if ~(  (ndims(y) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('y must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(z);
  if ~(  (ndims(z) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('z must be a vector, not e.g., meshgrid output');
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
    p = 3
    makeBanded = false;
  elseif (nargin == 7)
    makeBanded = false;
  elseif (nargin == 8)
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
  dz = z(2)-z(1);   Nz = length(z);
  ddx = [dx  dy  dz];
  ptL = [x(1) y(1) z(1)];

  if (Nx * Ny * Nz > 1e15)
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
  [Ibpt, Xgrid] = findGridInterpBasePt_test([xi yi zi], p, ptL, ddx);
  xw = LagrangeWeights1D_vec(Xgrid(:,1), xi, ddx(1), N);
  yw = LagrangeWeights1D_vec(Xgrid(:,2), yi, ddx(2), N);
  zw = LagrangeWeights1D_vec(Xgrid(:,3), zi, ddx(3), N);
  %toc

  %tic
  % this is a good order for memory access: ijk just counts up
  for k=1:N
    for i=1:N
      for j=1:N
        gi = (Ibpt(:,1) + i - 1);
        gj = (Ibpt(:,2) + j - 1);
        gk = (Ibpt(:,3) + k - 1);
        ijk = sub2ind([N,N,N], j, i, k);
        weights(:,ijk) = xw(:,i) .* yw(:,j) .* zw(:,k);

        % all these do the same, but last one is fastest.  Although sub2ind
        % presumably has safety checks...
        %ind = (gk-1)*(Nx*Ny) + (gi-1)*Ny + gj;
        %ind = sub2ind([Ny,Nx,Nz], gj, gi, gk);
        %ind = round((gk-1)*(Nx*Ny) + (gi-1)*(Ny) + gj-1 + 1);
        Ej(:,ijk) = (gk-1)*(Nx*Ny) + (gi-1)*Ny + gj;
      end
    end
  end
  %toc
  T1 = toc(T);
  %fprintf('done new Ei,Ej,weights, total time: %g\n', T1);


  % TODO: is there any advantage to keeping Ei as matrices?  Then each
  % column corresponds to the same point in the stencil...
  if ~makeListOutput
    tic
    E = sparse(Ei(:), Ej(:), weights(:), length(xi), Nx*Ny*Nz);
    T2 = toc;
    %fprintf('call to "sparse" time: %g\n', toc);
  end
  % Straightening them first doesn't make it faster
  %tic
  %Ei = Ei(:);
  %Ej = Ej(:);
  %weights = weights(:);
  %toc
  %tic
  %E = sparse(Ei, Ej, weights, length(xi), Nx*Ny*Nz);
  %toc

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
      Eout = E(:,setdiff(1:(Nx*Ny*Nz),band));
      if (nnz(Eout) > 0)
        nnz(Eout)
        warning('Lost some non-zero coefficients (from outside the innerband)');
      end
      E = Esparse;
      toc
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
  end
end
