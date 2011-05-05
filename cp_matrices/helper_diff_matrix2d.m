function L = helper_diff_matrix2d(x, y, band1, band2, weights, PTS, ndgrid)
%HELPER function: not intended for external use
%
% use_ndgrid: false for meshgrid, true for ndgrid
%
%   TODO: explain the roles of band1 and band2
%
%   TODO: this should work with dx != dy
%
%   TODO: issue with two bands: currently need extra padding in
%   wherever we call this from: should fix this


  % TODO why not just Nx = length(x)???
  %dx = x(2)-x(1);
  %dy = y(2)-y(1);
  %ddx = [dx  dy];
  %Nx = round( (x(end)-x(1)) / dx ) + 1;
  %Ny = round( (y(end)-y(1)) / dy ) + 1;
  %dim = length(ddx);
  Nx = length(x);
  Ny = length(y);

  StencilSize = length(weights);

  tic
  % Like a finite element code, find all the entries in 3 lists, then
  % insert them all at once while creating the sparse matrix
  Li = zeros((length(band1))*StencilSize, 1);
  Lj = zeros(size(Li));
  Ls = zeros(size(Li));
  Lc = 0;

  % main loop, a good candidate for parfor?
  for c = 1:length(band1)
    I = band1(c);

    if ndgrid
      % TODO: ndgrid ordering: not tested!
      [i,j] = ind2sub([Nx,Ny], I);
    else
      % meshgrid ordering
      [j,i] = ind2sub([Ny,Nx], I);
    end

    ii = i + PTS(:,1);
    jj = j + PTS(:,2);

    if ndgrid
      % TODO: ndgrid ordering: not tested!
      ind = round(sub2ind([Nx,Ny],ii,jj));
    else
      % funny ordering of y and x is b/c of meshgrid
      ind = round(sub2ind([Ny,Nx],jj,ii));
    end

    Lj( (Lc+1):(Lc+StencilSize) ) = ind;
    Li( (Lc+1):(Lc+StencilSize) ) = c*ones(size(ind));
    Ls( (Lc+1):(Lc+StencilSize) ) = weights;
    Lc = Lc + StencilSize;
  end

  if ( Lc ~= (length(band1)*StencilSize) )
    error('wrong number of elements');
  end

  L = sparse(Li, Lj, Ls, length(band1), Nx*Ny);
  L = L(:,band2);

  Ltime = toc;


