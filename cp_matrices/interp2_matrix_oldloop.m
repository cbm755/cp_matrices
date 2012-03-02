function E = interp2_matrix_oldloop(x, y, xi, yi, p, use_ndgrid)
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

  if (nargin < 6)
    use_ndgrid = false;
  end

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

    if use_ndgrid
      % ndgrid ordering
      % TODO: warning: maybe needs changes to buildInterpWeights?
      warning('not tested with ndgrid');
      ind = sub2ind([Nx,Ny], gii, gjj);
    else
      % meshgrid ordering
      ind = sub2ind([Ny,Nx], gjj, gii);
      %ind1 = round( (gii-1)*(Ny) + gjj-1 + 1 );
      %if (ind1 ~= ind)  error('fail')  end
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

  Etime = toc;
end
