function L = helper_diff_matrixnd(NN, band1, band2, weights, PTS)
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

  dim = length(NN);

  StencilSize = length(weights);

  tic
  % Like a finite element code, find all the entries in 3 lists, then
  % insert them all at once while creating the sparse matrix
  Li = repmat((1:length(band1))',1, StencilSize);
  Lj = zeros(size(Li));
  Ls = zeros(size(Li));

  ijk = myind2sub(NN, band1);  % todo: needs to match others?

  for c = 1:StencilSize
    for d=1:dim
      gi{d} = ijk(:,d) + PTS(c,d);
    end
    Lj(:,c) = sub2ind(NN, gi{:});
    Ls(:,c) = weights(c);
  end

  L = sparse(Li(:), Lj(:), Ls(:), length(band1), prod(NN));

  % remove columns not in band2 (the outerband)
  L = L(:, band2);

  Ltime = toc;
end


function A = myind2sub(siz, ndx)
% ndx is the linear index
  n = length(siz);
  k = [1 cumprod(siz(1:end-1))];
  A = zeros(length(ndx), n);
  for i = n:-1:1,
    vi = rem(ndx-1, k(i)) + 1;
    vj = (ndx - vi)/k(i) + 1;
    A(:,i) = vj;
    ndx = vi;
  end
end
