function [x,y,z] = paramHemisphere(n, R, cen)
%PARAMHEMISPHERE   A parameterization of a hemisphere
%   [x,y,z] = paramHemisphere(N, R) returns a mesh for a hemisphere of
%   radius R.  If R is omitted it defaults to 1.  The hemisphere is in
%   the (z>=0) half-space.  surf(x,y,z) can be used to make a plot.
%
%   [x,y,z] = paramHemisphere(N, R, CEN) returns a mesh centered
%   at 3-vector CEN.
%
%   Like Matlab's "SPHERE" but with a hemi:
%   http://www.youtube.com/watch?v=IyrcP5utXt4

  % defaults
  if (nargin < 2)
    R = 1;
  end
  if (nargin < 3)
    cen = [0,0,0];
  end

  % want to include the equator: easy way is to make sure we have an
  % even number of points when we call sphere()
  m = ceil(n/2);
  [xs,ys,zs] = sphere(2*m);
  x = xs( (m+1):(2*m+1), :);
  y = ys( (m+1):(2*m+1), :);
  z = zs( (m+1):(2*m+1), :);

  x = R*x + cen(1);
  y = R*y + cen(2);
  z = R*z + cen(3);
