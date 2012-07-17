function [E Ei Ej] = interp2_matrix_test(x, y, xi, yi, p)
%profile on
%INTERP2_MATRIX  Return a interpolation matrix
%   E = INTERP2_MATRIX(X,Y,XI,YI,DEGREEP)
%   Build a matrix which interpolates the grid data u onto the
%   points x and y using degree P barycentric Lagrange
%   interpolation.
%
%   E = INTERP2_MATRIX(X,Y,XI,YI,DEGREEP,true)
%   Same as above but assumes ndgrid rather than meshgrid ordering
%
%   Must be equispaced in each x,y (dx, dy can differ)
%   Currently does very little error checking of inputs

  % todo: in various places we assume xi is a vector
  %E = sparse(length(xi), length(u));
  
  % input checking
  [temp1,temp2] = size(x);
  if (temp1 ~= 1 && temp2 ~= 1)
    error('x must be a short vector, not meshgrid output');
  end
  [temp1,temp2] = size(y);
  if (temp1 ~= 1 && temp2 ~= 1)
    error('y must be a short vector, not meshgrid output');
  end
  [temp1,temp2] = size(xi);
  if (temp1 ~= 1 && temp2 ~= 1)
    error('yi must be a vector, not meshgrid output');
  end
  [temp1,temp2] = size(yi);
  if (temp1 ~= 1 && temp2 ~= 1)
    error('yi must be a vector, not meshgrid output');
  end
  xi = full(xi(:));
  yi = full(yi(:));
  
  N = p+1;
  % todo: assumes x, y are vectors
  ddx = [x(2)-x(1), y(2)-y(1)];
  dx = x(2)-x(1);
  dy = y(2)-y(1);
  Nx = round( (x(end)-x(1)) / dx ) + 1;
  Ny = round( (y(end)-y(1)) / dy ) + 1;

  if (Nx * Ny > 1e15)
    error('too big to use doubles as indicies: implement int64 indexing')
  end

  ptL = [x(1) y(1)];
  ptH = [x(end) y(end)];

  dim = length(ddx);
  ii = 0:N-1;
  jj = 0:N-1;
  [iig,jjg] = meshgrid(ii,jj);
  stencil = zeros(2,N^2);
  stencil(1,:) = iig(:)';
  stencil(2,:) = jjg(:)';
  
  EXTWIDTH = p+1;
  EXTSTENSZ = EXTWIDTH^dim;

  %tic
  
  Ei = repmat((1:length(xi))',1,EXTSTENSZ);
  weights = zeros(size(Ei));
  
  
  X = [xi yi];
  [out1, out2] = buildInterpWeights_test(X,ptL,ddx,p,stencil);
  xweights = out1{1};
  yweights = out1{2};
  gii = out2{1};
  gjj = out2{2};
  % meshgrid ordering
  Ej = sub2ind([Ny,Nx], gjj, gii); 
  
  
  for i = 1:N
      for j = 1:N
          weights(:,N*(i-1)+j) = xweights(:,i) .* yweights(:,j);
      end
  end
  
  
  Ei_vec = reshape(Ei,1,length(xi)*EXTSTENSZ);
  Ej_vec = reshape(Ej,1,length(xi)*EXTSTENSZ);
  Es = reshape(weights, 1, length(xi)*EXTSTENSZ);
  
  E = sparse(Ei_vec(:), Ej_vec(:), Es(:), length(xi), Nx*Ny);

  %Etime = toc
  %profile viewer
end
