function [x,y] = paramSimplecircle(n, R, cen)
%PARAMSEMICIRCLE   A parameterization of a semicircle

  % defaults
  if (nargin < 2)
    R = 1;
  end
  if (nargin < 3)
    cen = [0,0];
  end

  th = 0:2*pi/n:2*pi-2*pi/n;

  x = cos(th);
  y = sin(th);

  x = R*x + cen(1);
  y = R*y + cen(2);
