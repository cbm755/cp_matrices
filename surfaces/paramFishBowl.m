function [x,y,z] = paramHemisphereBand(n, R, cen)
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

  % default
  if (nargin < 2)
    R = 1;
  end
  if (nargin < 3)
    cen = [0,0,0];
  end

  % multiple by two is the easiest way to get the equator right.
  [xs,ys,zs] = sphere(2*n);
  zt = 0.4;

  ind = find(zs(:,1)<zt);

  x = xs(ind,:);
  y = ys(ind,:);
  z = zs(ind,:);

  x = R*x + cen(1);
  y = R*y + cen(2);
  z = R*z + cen(3);
  