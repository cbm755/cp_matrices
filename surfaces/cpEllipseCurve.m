function [cpx, cpy, dist, bdy] = cpEllipseCurve(x, y, lim, ab, cen)
%CPELLIPSE  Closest Point function for part of an ellipse.
%   [cpx, cpy, dist] = cpEllipseCurve(x, y, [lo up], [a b])
%      An ellipse centered at the origin with major axis 'a' and
%      minor axis 'b'.  'a' and 'b' default to 1.5 and 0.75.
%      The ellipse is cut-off in the x direction at 'lo' and 'up'.
%      TODO: currently measured in the 'theta' parameter of the
%      ellipse.
%      TODO: defaults
%   [cpx, cpy, dist] = cpEllipse(x, y, [lo up], [a b], cen)
%      The ellipse is centered at 'cen'.
%
%   Internally, calls cpEllipse().
%   TODO: only works for y >= 0
%
%   Note: returns unsigned distance

  % defaults
  if (nargin < 3) || isempty(ab)
    ab = [1.5 0.75];
  end
  a = ab(1);
  b = ab(2);
  if (nargin < 4) || isempty(lim)
    lim = [pi/10  4/5*pi];
  end
  if (nargin < 5) || isempty(cen)
    cen = [0 0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);

  % TODO: only works points which would be closest to the upper
  % half: points in the upper half plane
  if min(min(y)) < 0
    error('computation will be incorrect for bottom half plane')
  end

  leftx = a*cos(lim(2));
  rightx = a*cos(lim(1));
  lefty = b*sin(lim(2));
  righty = b*sin(lim(1));

  [cpx, cpy, dist] = cpEllipse(x, y, a, b);
  bdy = zeros(size(dist));
  dist = abs(dist);

  % fix x-axis points which point to the bottom
  I = cpy < 0;
  cpy(I) = abs(cpy(I));

  I = (cpx <= leftx);
  cpx(I) = leftx;
  cpy(I) = lefty;
  bdy(I) = 1;
  dist(I) = sqrt((x(I) - leftx).^2 + (y(I) - lefty).^2);

  I = (cpx >= rightx);
  cpx(I) = rightx;
  cpy(I) = righty;
  bdy(I) = 2;
  dist(I) = sqrt((x(I) - rightx).^2 + (y(I) - righty).^2);
  

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
