function E = interp3_matrix_oldloop(x, y, z, xi, yi, zi, p, use_ndgrid)
%INTERP3_MATRIX  Return a interpolation matrix
%   E = INTERP3_MATRIX(X,Y,Z,XI,YI,ZI,P)
%   Build a matrix which interpolates the grid data u onto the
%   points x and y using degree P barycentric Lagrange
%   interpolation.
%
%   E = INTERP2_MATRIX(X,Y,XI,YI,DEGREEP,true)
%   Same as above but assumes ndgrid rather than meshgrid ordering
%
%   Does very little error checking up the equispaced nature of x,y,z

  % todo: assumes xi is a vector, could relax this to return a
  % matrix sized based on the linear index (numels??).

  if (nargin < 8)
    use_ndgrid = false;
  end

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
  % 3D stencil size, p+1 points in each dimension
  StencilSize = (p+1)^dim;

  tic
  % old slow way is to allocate it first:
  %Eold = spalloc(usz+gsz, usz, StencilSize*(usz+gsz));

  % newer idea is like a finite element code, find all the entries
  % in 3 lists, then insert them all at once while creating the
  % sparse matrix
  Ei = zeros((length(xi))*StencilSize, 1);
  Ej = zeros(size(Ei));
  Es = zeros(size(Ei));
  Ec = 0;

  % good candidate for parfor?
  for i = 1:length(xi)
    X = [xi(i)  yi(i)  zi(i)];
    [weights,gii,gjj,gkk] = buildInterpWeights(X,ptL,ddx,p);

    if use_ndgrid
      % ndgrid ordering
      % TODO: warning: maybe needs changes to buildInterpWeights?
      warning('not tested with ndgrid');
      ind = sub2ind([Nx,Ny,Nz], gii, gjj, gkk);
    else
      % meshgrid ordering
      % (funny ordering of y and x)
      ind = sub2ind([Ny,Nx,Nz], gjj, gii, gkk);
      %ind = round((gkk-1)*(Nx*Ny) + (gii-1)*(Ny) + gjj-1 + 1);
    end

    %Eold(i,jj) = extWeights;
    % its faster to track all entries and put them all in at once
    Ej( (Ec+1):(Ec+StencilSize) ) = ind;
    Ei( (Ec+1):(Ec+StencilSize) ) = i*ones(size(ind));
    Es( (Ec+1):(Ec+StencilSize) ) = weights;
    Ec = Ec + StencilSize;
  end

  if ( Ec ~= (length(xi)*StencilSize) )
    error('wrong number of elements');
  end

  E = sparse(Ei, Ej, Es, length(xi), Nx*Ny*Nz);

  Etime = toc;
end

