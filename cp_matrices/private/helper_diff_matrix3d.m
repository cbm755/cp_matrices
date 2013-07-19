function L = helper_diff_matrix3d(x, y, z, band1, band2, weights, PTS, ndgrid)
%HELPER function: not intended for external use
%
%   TODO: explain the roles of band1 and band2
%
%   TODO: issue with two bands: currently need extra padding in
%   wherever we call this from
%
%   TODO: maybe pass an opt structure ndgrid, safety etc.

  % TODO: a global variable is probably not a good way to do this.
  global ICPM2009BANDINGCHECKS

  % less checks (e.g., negative indices)
  speedOverSafety = false;

  Nx = length(x);
  Ny = length(y);
  Nz = length(z);

  StencilSize = length(weights);

  tic
  % Like a finite element code, find all the entries in 3 lists, then
  % insert them all at once while creating the sparse matrix
  Li = repmat((1:length(band1))',1, StencilSize);
  Lj = zeros(size(Li));
  Ls = zeros(size(Li));

  if (ndgrid)
    [i,j,k] = ind2sub([Nx,Ny,Nz], band1);
  else
    [j,i,k] = ind2sub([Ny,Nx,Nz], band1);
  end

  for c = 1:StencilSize
    ii = i + PTS(c,1);
    jj = j + PTS(c,2);
    kk = k + PTS(c,3);
    Ls(:,c) = weights(c);
    if (ndgrid)
      if (speedOverSafety)
        Lj(:,c) = (kk-1)*(Nx*Ny) + (jj-1)*Nx + ii;
      else
        Lj(:,c) = sub2ind([Nx,Ny,Nz],ii,jj,kk);
      end
    else % meshgrid
      if (speedOverSafety)
        Lj(:,c) = (kk-1)*(Nx*Ny) + (ii-1)*Ny + jj;
      else
        Lj(:,c) = sub2ind([Ny,Nx,Nz],jj,ii,kk);
      end
    end

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
