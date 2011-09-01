function r = assertAlmostEqual(a,b,tol)
% ASSERTALMOSTEQUAL  Test for floating point equality 
%    r = assertAlmostEqual(a,b,tol)
%       returns true if abs(a-b) < tol
%       tol defaults to 1e-15 if left out.

  if (nargin < 3)
    tol = 100*eps;
  end

  r = all(abs(a-b) < tol);