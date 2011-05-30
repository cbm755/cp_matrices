function [cpx, cpy, dist] = cpCircle(x, y, R, xc, yc)
%CPCIRCLE  Closest Point function for a circle.
%   [cpx, cpy, dist] = cpCircle(x, y)
%      A unit circle centered at the origin
%   [cpx, cpy, dist] = cpCircle(x, y, R)
%      A circle of radius R centered at the origin
%   [cpx, cpy, dist] = cpCircle(x, y, R, xc, yc)
%      A circle of radius R centered at (xc,yc)

  % radius defaults to 1
  if (nargin < 3)
    R = 1;
  end
  if (nargin == 4)
    error('must specify both xc,yc');
  end
  % center shift, defaults to 0
  if (nargin < 5)
    xc = 0;
    yc = 0;
  end

  % shift to the origin
  x = x - xc;
  y = y - yc;

  [th, r] = cart2pol(x, y);
  [cpx, cpy] = pol2cart(th, R);

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 );

  % shift back
  cpx = cpx + xc;
  cpy = cpy + yc;


  % another approach:
  %r = sqrt(x.^2 + y.^2);
  %x2 = (r==0)*1.0 + (r~=0).*x;
  %r = (r==0)*1.0 + (r~=0).*r;
  %c = R ./ r;
  %cpx = c.*x2;
  %cpy = c.*y;
  %dist = sqrt( (cpx-x).^2 + (cpy-y).^2 );
