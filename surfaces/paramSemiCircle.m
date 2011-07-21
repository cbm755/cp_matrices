function [x,y] = paramSemiCircle(n, R, cen)
%PARAMHEMISPHERE   A parameterization of a hemisphere
%   [x,y,z] = paramHemisphere(N, R) returns a mesh for a hemisphere of
%   radius R.  If R is omitted it defaults to 1.  The hemisphere is in
%   the (z>=0) half-space.  surf(x,y,z) can be used to make a plot.
%
%   [x,y,z] = paramHemisphere(N, R, CEN) returns a mesh centered
%   at CEN.
%
%   N is twice as dense as you might think (mesh size is [N x 2*N]).
%
%   Like Matlab's "SPHERE" but with a hemi:
%   http://www.youtube.com/watch?v=IyrcP5utXt4

  % defaults
  if (nargin < 2)
    R = 1;
  end
  if (nargin < 3)
    cen = [0,0];
  end

  th = 0:pi/n:pi;
  
  x = cos(th);
  y = sin(th);
  
  x = R*x + cen(1);
  y = R*y + cen(2);
 