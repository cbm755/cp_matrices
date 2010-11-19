function E = interp2_matrix(x, y, xi, yi, p)
%INTERP2_MATRIX  Return a interpolation matrix
%   E = INTERP2_MATRIX(X,Y,XI,YI,DEGREEP)
%   Build a matrix which interpolates the grid data u onto the
%   points x and y using degree P barycentric Lagrange
%   interpolation.
%
%   Must be equispaced in each x,y (dx, dy can differ)
%   Does very little error checking of inputs

  % todo: here assumes xi is a vector
  %E = sparse(length(xi), length(u));

  % input checking
  [temp1,temp2] = size(x);
  if (temp1 ~= 1 & temp2 ~= 1)
    error('x must be a short vector, not meshgrid output');
  end
  [temp1,temp2] = size(y);
  if (temp1 ~= 1 & temp2 ~= 1)
    error('y must be a short vector, not meshgrid output');
  end
  [temp1,temp2] = size(xi);
  if (temp1 ~= 1 & temp2 ~= 1)
    error('yi must be a vector, not meshgrid output');
  end
  [temp1,temp2] = size(yi);
  if (temp1 ~= 1 & temp2 ~= 1)
    error('yi must be a vector, not meshgrid output');
  end

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
  EXTWIDTH = p+1;
  EXTSTENSZ = EXTWIDTH^dim;

  tic
  % old slow way
  %Eold = spalloc(usz+gsz, usz, PAR.EXTSTENSZ*(usz+gsz));

  Ei = zeros((length(xi))*EXTSTENSZ, 1);
  Ej = zeros(size(Ei));
  Es = zeros(size(Ei));
  Ec = 0;

  % good candidate for parfor?
  for i = 1:length(xi)
    X = [xi(i), yi(i)];
    [weights,gii,gjj] = buildInterpWeights(X,ptL,ddx,p);

    % can use sub2ind??
    % todo: Nx,Ny? order?
    % TODO: I hate this sort of thing, compare to CY's code
    ind1 = round( (gii-1)*(Ny) + gjj-1 + 1 );
    ind = round( sub2ind([Ny,Nx], gjj, gii) );
    if (ind1 ~= ind)
      error('fail')
    end

    %Eold(i,jj) = extWeights;
    % its faster to track all entries and put them all in at once
    Ej( (Ec+1):(Ec+EXTSTENSZ) ) = ind;
    Ei( (Ec+1):(Ec+EXTSTENSZ) ) = i*ones(size(ind));
    Es( (Ec+1):(Ec+EXTSTENSZ) ) = weights;
    Ec = Ec + EXTSTENSZ;
  end

  if ( Ec ~= (length(xi)*EXTSTENSZ) )
    error('wrong number of elements');
  end

  E = sparse(Ei, Ej, Es, length(xi), Nx*Ny);

  Etime = toc
end
