function [E,DiagwLinverse] = buildCPmatrixFromLaplacian3d(x, y, z, xi, yi, zi, p, band, use_ndgrid)
% buildCPmatrixFromLaplacian3d 
% Return a 3D pseudo-interpolation matrix with the help of discrete Laplacian.

% u(cp) = u(cp) + [\Delta u](cp) - [\Delta u](cp)
% we then approximate the first [\Delta u](cp) using finite differences,
% hoping to get some positive weights for off-diag entries, and replace the
% second [\Delta u](cp) by what we know from the orginal PDE.



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
  
  % relative coordinates w.r.t. the central block.
  hx = xi - Xgrid{1} - dx;
  hy = yi - Xgrid{2} - dx;
  hz = zi - Xgrid{3} - dx;
  
  % weights from diagonal entries of Laplacian.
  diagwL = 2 * ( 1./((dx+hx).*(2*dx-hx)) + 1./((dy+hy).*(2*dy-hy)) + ...
                1./((dz+hz).*(2*dz-hz)) );
            
  % weights from other 6 directions.
  xwL = cell(2,1); ywL = cell(2,1); zwL = cell(2,1);
  xwL{1} = 2/3/dx./(dx+hx);            % negative x-direction
  xwL{2} = 2/3/dx./(2*dx-hx);          % postive x-direction
  ywL{1} = 2/3/dy./(dy+hy);            % negative y-direction
  ywL{2} = 2/3/dy./(2*dy-hy);          % postive y-direction
  zwL{1} = 2/3/dz./(dz+hz);            % negative z-direction
  zwL{2} = 2/3/dz./(2*dz-hz);          % postive z-direction

  % compute the actual weights of the Laplacian at CP along the x-direction
  for cnt = 1:2
      % we want i = 1 & 4 correspondingly
      i = 3*cnt - 2;
      for k = 1:N
          for j = 1:N
              ijk = sub2ind([N,N,N], j, i, k);
              weights(:,ijk) = weights(:,ijk) + xwL{cnt}(:) .* yw(:,j) .* zw(:,k);
          end
      end
  end
  
  % compute the actual weights of the Laplacian at CP along the y-direction
  for cnt = 1:2
      % we want j = 1 & 4 correspondingly
      j = 3*cnt - 2;
      for k = 1:N
          for i = 1:N
              ijk = sub2ind([N,N,N], j, i, k);
              weights(:,ijk) = weights(:,ijk) + ywL{cnt}(:) .* xw(:,i) .* zw(:,k);
          end
      end
  end
  
  % compute the actual weights of the Laplacian at CP along the z-direction
  for cnt = 1:2
      % we want k = 1 & 4 correspondingly
      k = 3*cnt - 2;
      for i = 1:N
          for j = 1:N
              ijk = sub2ind([N,N,N], j, i, k);
              weights(:,ijk) = weights(:,ijk) + zwL{cnt}(:) .* xw(:,i) .* yw(:,j);
          end
      end
  end
  %toc

  
  
  %tic
  % compute the positions of the weights
  for k=1:N
    for i=1:N
      for j=1:N
        gi = (Ibpt{1} + i - 1);
        gj = (Ibpt{2} + j - 1);
        gk = (Ibpt{3} + k - 1);
        ijk = sub2ind([N,N,N], j, i, k);
        %weights(:,ijk) = xw(:,i) .* yw(:,j) .* zw(:,k);

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

  E = sparse(Ei(:), Ej(:), weights(:), length(xi), M);
  E = E(:,band);
  
  % divide by the diagonal entries of the Laplace matrix to actually get
  % the entries of E
  DiagwLinverse = spdiags([1./diagwL],0,length(xi),length(xi));
  %Epseudo = DiagwLinverse * E;
end
