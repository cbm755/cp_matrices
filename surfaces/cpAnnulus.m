function [cpx, cpy, dist, bdy] = cpAnnulus(x, y, R1, R2, cen)
%cpAnnulus  Closest Point function for an annulus
%   [cpx, cpy, dist, bdy] = cpAnnulus(x, y, R1, R2)
%   [cpx, cpy, dist, bdy] = cpAnnulus(x, y, R1, R2, cen)
%
%   R1 and R2 are the inner and outer radii.
%
%   The return value 'bdy' is 0 in the interior, 1 if the closest
%   point is on the interior boundary, and 2 if the closest point is
%   on the exterior boundary.
%
%   Examples:
%     >> [cpx, cpy, dist, bdy] = cpAnnulus(2, 0, 0.8, 1)
%     cpx = 1
%     cpy = 0
%     dist = 1
%     bdy = 2
%
%     >> [cpx, cpy, dist, bdy] = cpAnnulus(0.5, 0, 0.8, 1)
%     cpx = 0.80000
%     cpy = 0
%     dist = 0.30000
%     bdy = 1
%
%     >> [cpx, cpy, dist, bdy] = cpAnnulus(0.6, -0.6, 0.8, 1)
%     cpx = 0.60000
%     cpy = -0.60000
%     dist = 0
%     bdy = 0


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
