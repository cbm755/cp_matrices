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

  warning('only works for square case')
  warning('TODO: logical diagonal is not on diag of L in general!')
  % TODO: logical diagonals are the indices of band1 in band2?
  % extract the diagonal
  Ldiag = diag(diag(L));
  Ldiagpad = diag(diag(L),size(L,1),size(L,2));
  M = Ldiag + (L-Ldiagpad)*E;
