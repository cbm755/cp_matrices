function [x,y] = paramSimpleCircle(n, R, cen)
%PARAMSIMPLECIRCLE   A parameterization of a circle
%   [x, y] = paramSimpleCircle(n, R, cen)
%   [x, y, th] = paramSimpleCircle(n, R, cen)
%   'n' specifies the approximate number of points.
%   'R' and 'cen' are optional.

  % defaults
  if (nargin < 2)
    R = 1;
  end
  if (nargin < 3)
    cen = [0 0];
  end

  n = max(n, 3);

  %th = ( 0:(2*pi/n):(2*pi-2*pi/n) )';
  th = ( 0:(2*pi/n):(2*pi) )';

  x = cos(th);
  y = sin(th);

  x = R*x + cen(1);
  y = R*y + cen(2);
