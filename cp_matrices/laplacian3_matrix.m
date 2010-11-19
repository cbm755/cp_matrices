function L = laplacian3_matrix(x, y, z, order, band1, band2)
%LAPLACIAN3_MATRIX  Build a 3D discrete Laplacian
%   ORDER: 2 or 4 for 2nd or 4th-order
%   Does no error checking up the equispaced nature of x,y,z
% 
%   TODO: explain the roles of band1 and band2

  % input checking
  [temp1, temp2] = size(x);
  if ~(  (ndims(x) == 2) & (temp1 == 1 | temp2 == 1)  )
    error('x must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(y);
  if ~(  (ndims(y) == 2) & (temp1 == 1 | temp2 == 1)  )
    error('y must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(z);
  if ~(  (ndims(z) == 2) & (temp1 == 1 | temp2 == 1)  )
    error('z must be a vector, not e.g., meshgrid output');
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);
  ddx = [dx  dy  dz];
  Nx = round( (x(end)-x(1)) / dx ) + 1;
  Ny = round( (y(end)-y(1)) / dy ) + 1;
  Nz = round( (z(end)-z(1)) / dz ) + 1;

  ptL = [x(1) y(1) z(1)];
  ptH = [x(end) y(end) z(end)];

  dim = length(ddx);

  % TODO: code assumes dx=dy=dz?!
  if (order == 2)
    weights = [-6 1 1 1 1 1 1] / dx^2;
    PTS = [ 0   0   0; ...
            1   0   0; ...
           -1   0   0; ...
            0   1   0; ...
            0  -1   0; ...
            0   0   1; ...
            0   0  -1];
  elseif (order == 4)
    weights = [-15.0/2.0 ...
               (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
               (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
               (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
              ] / dx^2;
    PTS = [  0   0   0; ...
            -2   0   0; ...
            -1   0   0; ...
             1   0   0; ...
             2   0   0; ...
             0  -2   0; ...
             0  -1   0; ...
             0   1   0; ...
             0   2   0; ...
             0   0  -2; ...
             0   0  -1; ...
             0   0   1; ...
             0   0   2];
  else
    error(['order ' num2str(order) ' not implemented']);
  end
  StencilSize = length(weights);

  tic
  % Like a finite element code, find all the entries in 3 lists, then
  % insert them all at once while creating the sparse matrix
  Li = zeros((length(band1))*StencilSize, 1);
  Lj = zeros(size(Li));
  Ls = zeros(size(Li));
  Lc = 0;

  % good candidate for parfor?
  for c = 1:length(band1)
    I = band1(c);

    % meshgrid ordering
    [j,i,k] = ind2sub([Ny,Nx,Nz], I);

    %X = [x(i)  y(j)  z(k)];

    ii = i + PTS(:,1);
    jj = j + PTS(:,2);
    kk = k + PTS(:,3);

    % funny ordering of y and x is b/c of meshgrid
    ind = round(sub2ind([Ny,Nx,Nz],jj,ii,kk));

    Lj( (Lc+1):(Lc+StencilSize) ) = ind;
    Li( (Lc+1):(Lc+StencilSize) ) = c*ones(size(ind));
    Ls( (Lc+1):(Lc+StencilSize) ) = weights;
    Lc = Lc + StencilSize;
  end

  if ( Lc ~= (length(band1)*StencilSize) )
    error('wrong number of elements');
  end

  L = sparse(Li, Lj, Ls, length(band1), Nx*Ny*Nz);
  L = L(:,band2);

  Ltime = toc
end

