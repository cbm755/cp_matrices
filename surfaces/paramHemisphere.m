function [x,y,z] = paramHemisphere(n, varargin)
%HEMISPHERE  A parameterization of a hemisphere.
%   [x,y,z] = paramHemisphere(n, R) returns a mesh for a hemisphere of
%   radius R.  If R is omitted it defaults to 1.  The hemisphere is in
%   the (z>=0) half-space.  surf(x,y,z) can be used to make a plot.
%
%   [x,y,z] = paramHemisphere(n, R, xc,yc,zc) returns a mesh centered
%   at (xc,yc,zc).
% 
%   "n" is twice as dense as you might think (mesh size is [n x 2*n]).
%
%   Like Matlab's "SPHERE" but with a hemi:
%   http://www.youtube.com/watch?v=IyrcP5utXt4

  % multiple by two is the easiest way to get the equator right.
  [xs,ys,zs] = sphere(2*n);
  x = xs((n+1):2*n+1,:);
  y = ys((n+1):2*n+1,:);
  z = zs((n+1):2*n+1,:);

  if (nargin >= 2)
    R = varargin{1};
    if (nargin >= 3)
      xc = varargin{2};
      yc = varargin{3};
      zc = varargin{4};
    else
      xc = 0; yc = 0; zc = 0;
    end
    x = R*x + xc;
    y = R*y + yc;
    z = R*z + zc;
  end

