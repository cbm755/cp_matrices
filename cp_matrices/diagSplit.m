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

  % extract the diagonal
  Ldiag = diag(diag(L));
  M = Ldiag + (L-Ldiag)*E;
