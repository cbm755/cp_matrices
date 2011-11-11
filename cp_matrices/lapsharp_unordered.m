function M = lapsharp_unordered(L, E, R, delta)
%LAPSHARP_UNORDERED  Diagonal split matrix product of two matrices
%   M = LAPSHARP_UNORDERED(L, E, R)
%   Splits two matrices L and E as in Macdonald & Ruuth 2009.  L is
%   not square, and the "logical diagonal" is encoded in the
%   restriction operator R.
%
%   This generalized matrix product is for example, useful for
%   the implicit Closest Point Method, where is offers increased
%   stability.
%
%   M = LAPSHARP_UNORDERED(L, E, R, delta)
%   Laplacian sharp as in [Macdonald Brandman Ruuth 2011] (delta =
%   2*dim/epsilon^2).  Set delta to 2*dim/dx^2 (i.e., epsilon = dx)
%   to get the diagonal split operator.  Leave delta blank to
%   default to this behavior.


  if (nargin >= 4)
    %% Full lapsharp
    %delta = 2*dim/eps^2 = 2*dim/dx^2
    I = speye(size(L,1),size(E,2));
    M = L*E - delta*(I - R*E);
  else
    %% Diagonal splitting
    % a particular case of the above (implemented with less
    % rounding error)

    %logical diagonal is not on diag of L in general!  (See the
    %implementation in that case in lapsharp.m)

    [i,j,r] = find(R);  % indices of diagonal elements
    Ldiagpad = R .* L;
    Ldiag = Ldiagpad(i,j);  % a square matrix
    M = Ldiag + (L-Ldiagpad)*E;

    % this gives error for large matrix
    %dd = find(R);       % linear index of diagonal
    %Ldiag = diag(L(dd));
    %Ldiagpad = sparse(i,j,L(dd),size(L,1),size(L,2));
    % No, its not this:
    % Ldiagpad = spdiags(L(dd), 0, size(L,1), size(L,2));
  end