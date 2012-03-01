function L = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, ndgrid)
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

  if (ndgrid)
    error('ndgrid not implemented in new vector code');
    % TODO: fix this, or remove this option
  end

  Nx = length(x);
  Ny = length(y);
  Nz = length(z);

  StencilSize = length(weights);

  tic
  % Like a finite element code, find all the entries in 3 lists, then
  % insert them all at once while creating the sparse matrix
  %  Li = zeros((length(band1))*StencilSize, 1);
  %Li = zeros((length(band1)),StencilSize);
  Li = repmat((1:length(band1))',1,StencilSize);
  Lj = zeros(size(Li));
  Ls = zeros(size(Li));
  Lc = 0;

  [j,i,k] = ind2sub([Ny,Nx,Nz], band1);
  for c = 1:StencilSize
    ii = i + PTS(c,1);
    jj = j + PTS(c,2);
    kk = k + PTS(c,3);
    % TODO: this is faster but no satefy checks, perhaps we need a
    % "safety" option...
    %Lj(:,c) = (kk-1)*(Nx*Ny) + (ii-1)*Ny + jj;
    Lj(:,c) = sub2ind([Ny,Nx,Nz],jj,ii,kk);
    Ls(:,c) = weights(c);
  end

  L = sparse(Li(:), Lj(:), Ls(:), length(band1), Nx*Ny*Nz);

  % TODO: these sorts of checks could move to the ops and bands replacement
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
