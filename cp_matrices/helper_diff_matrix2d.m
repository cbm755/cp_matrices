function L = helper_diff_matrix2d(x, y, band1, band2, weights, PTS, ndgrid)
%HELPER function: not intended for external use
%
%   TODO: explain the roles of band1 and band2
%
%   TODO: issue with two bands: currently need extra padding in
%   wherever we call this from
%
%   TODO: maybe pass an opt structure ndgrid, safety etc.

  % TODO: global variable may not be the best way to do this
  global ICPM2009BANDINGCHECKS

  % less checks (e.g., negative indices)
  speedOverSafety = false;

  Nx = length(x);
  Ny = length(y);

  StencilSize = length(weights);

  tic
  % Like a finite element code, find all the entries in 3 lists, then
  % insert them all at once while creating the sparse matrix
  Li = repmat((1:length(band1))', 1, StencilSize);
  Lj = zeros(size(Li));
  Ls = zeros(size(Li));

  if (ndgrid)
    [i,j] = ind2sub([Nx,Ny], band1);
  else
    [j,i] = ind2sub([Ny,Nx], band1);
  end

  for c = 1:StencilSize
    ii = i + PTS(c,1);
    jj = j + PTS(c,2);
    Ls(:,c) = weights(c);
    if (ndgrid)
      if (speedOverSafety)
        Lj(:,c) = (jj-1)*Nx + ii;
      else
        Lj(:,c) = sub2ind([Nx,Ny],ii,jj);
      end
    else % meshgrid
      if (speedOverSafety)
        Lj(:,c) = (ii-1)*Ny + jj;
      else
        Lj(:,c) = sub2ind([Ny,Nx],jj,ii);
      end
    end
  end
  L = sparse(Li(:), Lj(:), Ls(:), length(band1), Nx*Ny);

  % TODO: these sorts of checks could move to the ops and bands replacement
  % If we're using careful banding a la iCPM2009 then as a sanity
  % check all of the columns outside of band2 should be zero.
  if (~isempty(ICPM2009BANDINGCHECKS)) && (ICPM2009BANDINGCHECKS)
    Lout = L(:, setdiff(1:(Nx*Ny),band2));
    if (nnz(Lout) > 0)
      nnz(Lout)
      error('Lost some non-zero coefficients (from outside the outerband)');
    end
  end

  % remove columns not in band2 (the outerband)
  L = L(:, band2);

  Ltime = toc;
end
