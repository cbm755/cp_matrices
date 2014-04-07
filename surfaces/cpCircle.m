function [cpx, cpy, sdist] = cpCircle(x, y, R, cen)
%CPCIRCLE  Closest Point function for a circle.
%   [cpx, cpy, sdist] = cpCircle(x, y)
%      A unit circle centered at the origin.
%   [cpx, cpy, sdist] = cpCircle(x, y, R)
%      A circle of radius R centered at the origin.
%   [cpx, cpy, sdist] = cpCircle(x, y, R, CEN)
%      A circle of radius R centered at CEN = [xc,yc].
%
%   Note: returns signed distance (with negative inside).


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

  sdist = r - R;

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);


  % another approach:
  %r = sqrt(x.^2 + y.^2);
  %x2 = (r==0)*1.0 + (r~=0).*x;
  %r = (r==0)*1.0 + (r~=0).*r;
  %c = R ./ r;
  %cpx = c.*x2;
  %cpy = c.*y;
  %dist = sqrt( (cpx-x).^2 + (cpy-y).^2 );
