function [x, y] = paramRoundedSquare(n, radii, cen)
%PARAMROUNDEDSQUARE  Parameterization of a square with rounded corners
%   [x, y] = paramRoundedSquare(n, radii, center)
%   'n' is the approximate number of points.

  % defaults
  if (nargin < 2)
    radii = 0.25;
  end
  if (nargin < 3)
    cen = [0 0];
  end

  if (length(radii) == 1)
    radii = radii*[1 1 1 1];
  end

  paramfs = {};
  paramfs{1} = @(n) paramArc(n, radii(1), [ 1-radii(1)  1-radii(1)], 0,      pi/2);
  paramfs{3} = @(n) paramArc(n, radii(2), [-1+radii(2)  1-radii(2)], pi/2,   pi);
  paramfs{5} = @(n) paramArc(n, radii(3), [-1+radii(3) -1+radii(3)], pi,     3*pi/2);
  paramfs{7} = @(n) paramArc(n, radii(4), [ 1-radii(4) -1+radii(4)], 3*pi/2, 0);

  Pts = {};
  Pts{1} = [ 1-radii(1)   1];
  Pts{2} = [-1+radii(2)   1];
  Pts{3} = [-1   1-radii(2)];
  Pts{4} = [-1  -1+radii(3)];
  Pts{5} = [-1+radii(3)  -1];
  Pts{6} = [ 1-radii(4)  -1];
  Pts{7} = [ 1  -1+radii(4)];
  Pts{8} = [ 1   1-radii(1)];
  paramfs{2} = @(n) paramLineSegment2d(n, Pts{1}, Pts{2});
  paramfs{4} = @(n) paramLineSegment2d(n, Pts{3}, Pts{4});
  paramfs{6} = @(n) paramLineSegment2d(n, Pts{5}, Pts{6});
  paramfs{8} = @(n) paramLineSegment2d(n, Pts{7}, Pts{8});

  % how many points per segment?  first compute the arclenths
  arclens = [radii(1)*pi/2  ...
             (2 - radii(1) - radii(2)) ...
             radii(2)*pi/2  ...
             (2 - radii(2) - radii(3))  ...
             radii(3)*pi/2  ...
             (2 - radii(2) - radii(3))  ...
             radii(4)*pi/2 ...
             (2 - radii(4) - radii(1))];
  arclens = arclens / sum(arclens);
  m = ceil(arclens * n);

  [x,y] = paramfs{1}(m(1));
  for j=2:length(paramfs)
    [x2,y2] = paramfs{j}(m(j));
    x = [x; x2(2:end)];
    y = [y; y2(2:end)];
  end

  x = x + cen(1);
  y = y + cen(2);
