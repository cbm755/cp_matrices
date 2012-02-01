function M = lapsharp(L, E, delta)
%LAPSHARP  Diagonal split matrix product of two matrices
%   M = LAPSHARP(L, E)
%   Splits two matrices L and E as in Macdonald & Ruuth 2009.
%   This generalized matrix product is for example, useful for
%   the implicit Closest Point Method, where is offers increased
%   stability.
%
%   See also LAPSHARP_UNORDERED

  if (nargin >= 3)
    %% Full lapsharp
    %delta = 2*dim/eps^2 = 2*dim/dx^2
    warning('not really tested');
    I = speye(size(L,1),size(E,2));
    M = L*E - delta*(I - E);

  else
    %% Diagonal splitting (default)
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
  end
