function [x,y,z] = paramSphere(n, varargin)
%PARAMSPHERE  A parameterization of a sphere.
%   [x,y,z] = paramSphere(n, R) returns a mesh for a sphere of radius
%   R.  If R is omitted it defaults to 1.  surf(x,y,z) can be used to
%   make a plot.
%
%   [x,y,z] = paramSphere(n, R, xc,yc,zc) returns a mesh for a sphere
%   of radius R centered at (xc,yc,zc).
%
%   A wrapper of Matlab's "SPHERE".

  [x,y,z] = sphere(n);

  if (nargin >= 2)
    R = varargin{1}
    if (nargin >= 3)
      xc = varargin{2}
      yc = varargin{3}
      zc = varargin{4}
    else
      xc = 0; yc = 0; zc = 0;
    end
    x = R*x + xc;
    y = R*y + yc;
    z = R*z + zc;
  end

