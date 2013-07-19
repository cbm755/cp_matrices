function L = helper_diff_matrix3d_oldloop(x, y, z, band1, band2, weights, PTS, ndgrid)
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

  % TODO: a global variable is probably not a good way to do this.
  global ICPM2009BANDINGCHECKS

  Nx = length(x);
  Ny = length(y);
  Nz = length(z);

  StencilSize = length(weights);

  tic
  % Like a finite element code, find all the entries in 3 lists, then
  % insert them all at once while creating the sparse matrix
  Li = zeros((length(band1))*StencilSize, 1);
  Lj = zeros(size(Li));
  Ls = zeros(size(Li));
  Lc = 0;

  % main loop, good candidate for parfor?
  for c = 1:length(band1)
    I = band1(c);

    if ndgrid
      % TODO: ndgrid ordering, not tested
      [i,j,k] = ind2sub([Nx,Ny,Nz], I);
    else
      % meshgrid ordering
      [j,i,k] = ind2sub([Ny,Nx,Nz], I);
    end

    ii = i + PTS(:,1);
    jj = j + PTS(:,2);
    kk = k + PTS(:,3);

    if ndgrid
      % TODO: ndgrid ordering: not tested!
      ind = sub2ind([Nx,Ny,Nz],ii,jj,kk);
    else
      % meshgrid ordering
      % used to round() here, could also think about using some integer type
      ind = sub2ind([Ny,Nx,Nz],jj,ii,kk);
    end

    Lj( (Lc+1):(Lc+StencilSize) ) = ind;
    Li( (Lc+1):(Lc+StencilSize) ) = c*ones(size(ind));
    Ls( (Lc+1):(Lc+StencilSize) ) = weights;
    Lc = Lc + StencilSize;
  end

  if ( Lc ~= (length(band1)*StencilSize) )
    error('wrong number of elements');
  end

  L = sparse(Li, Lj, Ls, length(band1), Nx*Ny*Nz);

  % If we're using careful banding a la iCPM2009 then as a sanity
  % check all of the columns outside of band2 should be zero.
  if (~isempty(ICPM2009BANDINGCHECKS)) && (ICPM2009BANDINGCHECKS)
    Lout = L(:, setdiff(1:(Nx*Ny*Nz),band2));
    if (nnz(Lout) > 0)
      nnz(Lout)
      error('Lost some non-zero coefficients (from outside the outerband)');
    end
  end

  % remove columns not in band2 (the outerband)
  L = L(:, band2);

  Ltime = toc;
end
