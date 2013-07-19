function L = helper_diff_matrixnd(NN, band1, band2, weights, PTS, invbandmap)
%HELPER function: not intended for external use
%
%   TODO: expose the speed versus memory tradeoff?  Default should
%   be less memory because it can be quite dramatic in say 5D.

  SaveMemOverSpeed = true;

  dim = length(NN);

  StencilSize = length(weights);

  % Like a finite element code, find all the entries in 3 lists, then
  % insert them all at once while creating the sparse matrix
  Li = repmat((1:length(band1))',1, StencilSize);
  Lj = zeros(size(Li));
  Ls = zeros(size(Li));


  [ijk{1:dim}] = ind2sub(NN, band1);


  if SaveMemOverSpeed
    %% efficient use of memory, but slightly slower

    if isempty(invbandmap)
      invbandmap = make_invbandmap(prod(NN), band2);
    end

    for c = 1:StencilSize
      for d=1:dim
        gi{d} = ijk{d} + PTS(c,d);
        I = find((gi{d} <= 0) | (gi{d} > NN(d)));
        if ~isempty(I)
          error('currently requires some padding between band and box');
        end
      end
      % sub2ind here needs to find all the stencil points in the
      % box of size NN.
      Lj(:,c) = invbandmap(sub2ind(NN, gi{:}));
      Ls(:,c) = weights(c);
    end

    I = (Lj ~= 0);
    L = sparse(Li(I), Lj(I), Ls(I), length(band1), length(band2));
    % For Laplacian, discarded values should not include diagonals
    % so no negative values (e.g., -2 in 1,-2,1).  To do this in
    % general might be harder.
    %assert(min(Ls(~I)) > 0);

  else
    %% faster but possibly memory wasteful implementation
    for c = 1:StencilSize
      for d=1:dim
        gi{d} = ijk{d} + PTS(c,d);
      end
      Lj(:,c) = sub2ind(NN, gi{:});
      Ls(:,c) = weights(c);
    end
    L = sparse(Li(:), Lj(:), Ls(:), length(band1), prod(NN));
    % remove columns not in band2 (the outerband)
    L = L(:, band2);
  end

