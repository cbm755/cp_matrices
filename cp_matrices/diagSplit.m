function M = diagSplit(L, E)
%DIAGSPLIT  Diagonal split matrix product of two matrices
%   M = DIAGSPLIT(L, E)
%   Splits two matrices L and E as in Macdonald & Ruuth 2009.
%   This generalized matrix product is for example, useful for
%   the implicit Closest Point Method, where is offers increased
%   stability.
%
%   TODO: Is this the most general implementation?  Does it work for
%   nonsquare L and E?

  if (size(L,1) ~= size(L,2))
    warning('*** only works for square case ***')
  end
  %warning('TODO: logical diagonal is not on diag of L in general!')
  % TODO: logical diagonals are the indices of band1 in band2?
  % extract the diagonal

  % surprisingly, diag here gives a sparse matrix
  Ldiag = diag(diag(L));
  Ldiagpad = spdiags(diag(L), 0, size(L,1), size(L,2));
  M = Ldiag + (L-Ldiagpad)*E;
