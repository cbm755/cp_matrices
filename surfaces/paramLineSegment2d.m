function [x, y] = paramLineSegment2d(n, p, q)
%PARAMLINESEGMENT2D  Parameterization of a line segment in 2D
%   [x, y] = paraLineSegment2d(n, p, q)
%   'p' and 'q' are two points in 2D.  'n' is the number of
%   output points.

  if (nargin == 1)
    p = [0 0];
    q = [1 0];
  end

  x = linspace(p(1), q(1), max(n,2))';
  y = linspace(p(2), q(2), max(n,2))';
