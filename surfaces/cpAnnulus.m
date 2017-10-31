function [cpx, cpy, dist, bdy] = cpAnnulus(x, y, R1, R2, cen)
%CPDISC  Closest Point function for a annulus
%   [cpx, cpy, dist, bdy] = cpAnnulus(x, y, R1, R2)


  % defaults
  if (nargin < 3)
    R1 = 1/2;
  end
  if (nargin < 4)
    R2 = 1;
  end
  if (nargin < 5)
    cen = [0 0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);

  [th, r] = cart2pol(x, y);
  [cpx1, cpy1] = pol2cart(th, R1);
  [cpx2, cpy2] = pol2cart(th, R2);

  dist1 = r - R1;  % negative inside R1
  dist2 = r - R2;  % negative inside R2

  bdy1 = dist1 < 0;
  bdy2 = dist2 > 0;
  cpx = cpx1;
  cpy = cpy1;

  cpx(bdy2) = cpx2(bdy2);
  cpy(bdy2) = cpy2(bdy2);

  I = (dist1 > 0) & (dist2 < 0);
  cpx(I) = x(I);
  cpy(I) = y(I);
  dist = dist1;
  dist(I) = 0;

  bdy = bdy1*1 + bdy2*2;

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
