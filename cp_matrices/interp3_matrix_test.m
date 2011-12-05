function E = interp3_matrix_test(x, y, z, xi, yi, zi, p)
%profile on

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
  [temp1, temp2] = size(xi);
  if ~(  (ndims(xi) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('xi must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(yi);
  if ~(  (ndims(yi) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('yi must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(zi);
  if ~(  (ndims(zi) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('zi must be a vector, not e.g., meshgrid output');
  end
  
  xi = full(xi(:));
  yi = full(yi(:));
  zi = full(zi(:));
  
  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);
  ddx = [dx  dy  dz];
  Nx = round( (x(end)-x(1)) / dx ) + 1;
  Ny = round( (y(end)-y(1)) / dy ) + 1;
  Nz = round( (z(end)-z(1)) / dz ) + 1;

  if (Nx * Ny * Nz > 1e15)
    error('too big to use doubles as indicies: implement int64 indexing')
  end

  ptL = [x(1) y(1) z(1)];
  ptH = [x(end) y(end) z(end)];

  dim = length(ddx);
  N = p+1;
  EXTSTENSZ = N^dim;

  if (Nx * Ny > 1e15)
    error('too big to use doubles as indicies: implement int64 indexing')
  end
  
  ii = 0:N-1;
  jj = 0:N-1;
  kk = 0:N-1;
  [iig,jjg,kkg] = meshgrid(ii,jj,kk);
  stencil = zeros(3,N^dim);
  stencil(1,:) = iig(:)';
  stencil(2,:) = jjg(:)';
  stencil(3,:) = kkg(:)';
  
  
  tic
  
  Ei = repmat((1:length(xi))',1,EXTSTENSZ);
  weights = zeros(size(Ei));
    
  X = [xi yi zi];
  [out1, out2] = buildInterpWeights_test(X,ptL,ddx,p,stencil);
  xweights = out1{1};
  yweights = out1{2};
  zweights = out1{3};
  gii = out2{1};
  gjj = out2{2};
  gkk = out2{3};
  clear out1 out2

  % meshgrid ordering
  Ej = sub2ind([Ny,Nx,Nz], gjj, gii, gkk);
  clear gii gjj gkk
  
  for i = 1:N
      for j = 1:N
          for k = 1:N
              weights(:,N^2*(k-1)+N*(i-1)+j) = xweights(:,i) .* yweights(:,j) .* zweights(:,k);
          end
      end
  end
  
  %Ei = reshape(Ei,1,length(xi)*EXTSTENSZ);
  %Ej = reshape(Ej,1,length(xi)*EXTSTENSZ);
  %Es = reshape(weights, 1, length(xi)*EXTSTENSZ);
  E = sparse(Ei(:), Ej(:), weights(:), length(xi), Nx*Ny*Nz);
  clear Ei Ej weights xweights yweights zweights

  Etime = toc
  %profile viewer
end
