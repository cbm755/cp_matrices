function [cpxx, cpyy, sdist] = cpSquare(xx, yy, cen)
%CPSQUARE  Closest Point function for a square
%   [cpx, cpy, sdist] = cpSquare(x, y)
%      A square with side length 2 centered centered at
%      the origin.
%   [cpx, cpy, sdist] = cpSquare(x, y, xc, yc)
%      A square with side length 2 centered centered at
%      the point (xc,yc)
%
%   Note: returns signed distance (with negative inside).
%
%   Not vectorized, uses a loop internally.

  % defaults
  if (nargin < 3)
    cen = [0, 0];
  end

  % shift to the origin
  xx = xx - cen(1);
  yy = yy - cen(2);

  cpxx = zeros(size(xx));
  cpyy = zeros(size(yy));

  for i=1:length(xx(:))
    x = xx(i);
    y = yy(i);

    if (y >= 1)
      cpy = 1;
    elseif (y <= -1)
      cpy = -1;
    elseif (y >= abs(x))
      cpy = 1;
    elseif (y <= -abs(x))
      cpy = -1;
    else
      cpy = y;
    end

    if (x >= 1)
      cpx = 1;
    elseif (x <= -1)
      cpx = -1;
    elseif (x > abs(y))
      cpx = 1;
    elseif (x < -abs(y))
      cpx = -1;
    else
      cpx = x;
    end

    cpxx(i) = cpx;
    cpyy(i) = cpy;
  end

  dist = sqrt( (xx-cpxx).^2 + (yy-cpyy).^2 );
  sdist = max( abs(xx), abs(yy) ) - 1;
  sdist = (sdist>0) .* max(dist, sdist)  +   (sdist <= 0) .* sdist;

  % shift back
  cpxx = cpxx + cen(1);
  cpyy = cpyy + cen(2);

