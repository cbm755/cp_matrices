function [varargout] = paramLineSegment(n, p, q)
%PARAMLINESEGMENT  Parameterization of a line segment in 2D/3D
%   [x,y] = paraLineSegment(n, p, q)
%   [x,y,z] = paraLineSegment(n, p, q)
%   'p' and 'q' are two points in 2D/3D.  'n' is the number of
%   output points.

  dim = length(p);

  varargout = {};
  for j=1:dim
    varargout{j} = linspace(p(j), q(j), max(n,2));
  end
