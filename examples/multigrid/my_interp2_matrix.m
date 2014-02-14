function [E,Ej,Es] = my_interp2_matrix(x, y, xi, yi, p, band, use_ndgrid)

  %% linear interpolation using cubic stencil. 
  
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
    use_ndgrid = false;
  elseif (nargin == 5)
    makeBanded = false;
    use_ndgrid = false;
  elseif (nargin == 6)
    makeBanded = true;
    use_ndgrid = false;
  elseif (nargin == 7)
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
  N = 2;
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
  [Ibpt, Xgrid] = findGridInterpBasePt_vec([xi yi], 3, ptL, ddx);
  xw = LagrangeWeights1D_vec(Xgrid(:,1), xi, 3*ddx(1), N);
  yw = LagrangeWeights1D_vec(Xgrid(:,2), yi, 3*ddx(2), N);
  %toc

  %tic
  % this is a good order for memory access: ijk just counts up
  for i=1:N
    for j=1:N
      gi = (Ibpt(:,1) + i - 1);
      gj = (Ibpt(:,2) + j - 1);
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
