function [cpx, cpy, dist, bdy] = cpDisc(x, y, R, cen)
%CPDISC  Closest Point function for a disc.
%   [cpx, cpy, dist, bdy] = cpDisc(x, y)
%      A unit disc centered at the origin.
%

  % defaults
  if (nargin < 3)
    R = 1;
  end
  if (nargin < 4)
    cen = [0 0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);

  [th, r] = cart2pol(x, y);
  [cpx, cpy] = pol2cart(th, R);

  dist = r - R;

  I = dist < 0;
  cpx(I) = x(I);
  cpy(I) = y(I);
  dist(I) = 0;

  bdy = ~I;

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
