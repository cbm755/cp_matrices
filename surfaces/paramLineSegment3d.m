function [x, y, z] = paramLineSegment3d(n, p, q)
%PARAMLINESEGMENT3D  Parameterization of a line segment in 3D
%   [x, y, z] = paraLineSegment3d(n, p, q)
%   'p' and 'q' are two points in 3D.  'n' is the number of
%   output points.

  if (nargin == 1)
    p = [0 0 0];
    q = [1 0 0];
  end

  x = linspace(p(1), q(1), max(n,2))';
  y = linspace(p(2), q(2), max(n,2))';
  z = linspace(p(3), q(3), max(n,2))';
