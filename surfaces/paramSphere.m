function [x,y,z] = paramSphere(n, R, cen)
%PARAMSPHERE  A parameterization of a sphere
%   [x,y,z] = paramSphere(N, R) returns a mesh for a sphere of radius
%   R.  If R is omitted it defaults to 1.  surf(x,y,z) can be used to
%   make a plot.
%
%   [x,y,z] = paramSphere(N, R, CEN) returns a mesh centered at CEN.
%
%   A wrapper of Matlab's "SPHERE".

  % defaults
  if (nargin < 2)
    R = 1;
  end
  if (nargin < 3)
    cen = [0,0,0];
  end

  [x,y,z] = sphere(n);

  x = R*x + cen(1);
  y = R*y + cen(2);
  z = R*z + cen(3);
