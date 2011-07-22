function [x,y] = paramEllipse(n, aa,bb, cen)
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
    aa = 1;
    bb = 1;
  end
  if (nargin < 3)
    cen = [0,0];
  end

  
  th = 0:2*pi/n:2*pi-2*pi/n;
  x = aa*cos(th);
  y = bb*sin(th);
