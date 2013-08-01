function [cpx, cpy, dist] = cpSimpleCircle(x, y, R)
%CPSIMPLECIRCLE  Closest Point function for a circle.
%   [cpx, cpy, dist] = cpCircle(x, y)
%      A unit circle centered at the origin
%   [cpx, cpy, dist] = cpCircle(x, y, R)
%      A circle of radius R centered at the origin

  % default value for R
  if (nargin < 3)
    R = 1;
  end

  % convert to polar and then back with a new radius
  [th, r] = cart2pol(x, y);
  [cpx, cpy] = pol2cart(th, R);

  % distance from (x,y) to (cpx,cpy)
  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 );
