function [E,GAMMA] = my_interp3_matrix(x, y, z, xi, yi, zi, p, band, use_ndgrid)
% hope to construct a 4th-order extension matrix without using corner
% points of the 4*4*4 stencil.
% also return the corresponding GAMMA matrix so that E1*L - GAMMA*(I-E) is
% an m-matrix.

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
  xw3 = LagrangeWeights1D_vec(Xgrid{1}, xi, ddx(1), N);
  yw3 = LagrangeWeights1D_vec(Xgrid{2}, yi, ddx(2), N);
  zw3 = LagrangeWeights1D_vec(Xgrid{3}, zi, ddx(3), N);
  
  xw1 = LagrangeWeights1D_vec(Xgrid{1}+dx, xi, ddx(1), 2);
  yw1 = LagrangeWeights1D_vec(Xgrid{2}+dy, yi, ddx(2), 2);
  zw1 = LagrangeWeights1D_vec(Xgrid{3}+dz, zi, ddx(3), 2);
   
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
        if i == 2 || i == 3
            weights(:,ijk) = xw1(:,i-1) .* yw3(:,j) .* zw3(:,k);
        end
        
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
  
  for j = 2:3
      for k = 2:3
          ijk = sub2ind([N,N,N], j, 1 ,k);
          weights(:,ijk) = weights(:,ijk) + yw1(:,j-1) .* zw1(:,k-1) .* xw3(:,1);
      end
  end
  for j = 2:3
      for k = 2:3
          ijk = sub2ind([N,N,N], j, 2 ,k);
          weights(:,ijk) = weights(:,ijk) - 2*yw1(:,j-1) .* zw1(:,k-1) .* xw3(:,1);
      end
  end
  for j = 2:3
      for k = 2:3
          ijk = sub2ind([N,N,N], j, 3 ,k);
          weights(:,ijk) = weights(:,ijk) + yw1(:,j-1) .* zw1(:,k-1) .* xw3(:,1);
      end
  end
  
  for j = 2:3
      for k = 2:3
          ijk = sub2ind([N,N,N], j, 2 ,k);
          weights(:,ijk) = weights(:,ijk) + yw1(:,j-1) .* zw1(:,k-1) .* xw3(:,4);
      end
  end
  for j = 2:3
      for k = 2:3
          ijk = sub2ind([N,N,N], j, 3 ,k);
          weights(:,ijk) = weights(:,ijk) - 2*yw1(:,j-1) .* zw1(:,k-1) .* xw3(:,4);
      end
  end
  for j = 2:3
      for k = 2:3
          ijk = sub2ind([N,N,N], j, 4 ,k);
          weights(:,ijk) = weights(:,ijk) + yw1(:,j-1) .* zw1(:,k-1) .* xw3(:,4);
      end
  end
  
  
  %toc
  T1 = cputime() - T1;
  %fprintf('done new Ei,Ej,weights, total time: %g\n', T1);


  % TODO: is there any advantage to keeping Ei as matrices?  Then each
  % column corresponds to the same point in the stencil...
 
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


  b = yi - Xgrid{2} - dx;
  c = zi - Xgrid{3} - dx;
  gamma = 24*dx^2 ./ ( (dx+b).*(2*dx-b).*(dx+c).*(2*dx-c) );
  GAMMA = spdiags(gamma,[0],length(xi),length(xi));
  
 
  end
