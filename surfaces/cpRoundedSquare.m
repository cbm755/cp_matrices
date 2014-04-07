function [cpx, cpy, sdist] = cpRoundedSquare(x, y, radii, cen)
%CPROUNDEDSQUARE  Closest Point function for a square with rounded corners
%   [cpx, cpy, dist] = cpRoundedSquare(x, y, radius)
%      A square with side length 2 centered centered at
%      the origin, each corner is rounded with 'radius'
%   [cpx, cpy, dist] = cpRoundedSquare(x, y, radii)
%      A square with side length 2 centered centered at
%      the origin, the corner of the jth quadrant is rounded with
%      radius 'radii(j)'.
%   [cpx, cpy, dist] = cpRoundedSquare(x, y, radii, [xc yc])
%      As above, centered centered at the point (xc,yc)
%
%   Note: returns signed distance (with negative inside).


  % defaults
  if (nargin < 3)
    radii = 0.25;
  end
  if (nargin < 4)
    cen = [0 0];
  end

  if (length(radii) == 1)
    radii = radii*[1 1 1 1];
  end

  x = x - cen(1);
  y = y - cen(2);

  cpfs = {};
  cpfs{1} = @(x,y) cpArc(x,y, radii(1), [ 1-radii(1)  1-radii(1)], 0,      pi/2);
  cpfs{3} = @(x,y) cpArc(x,y, radii(2), [-1+radii(2)  1-radii(2)], pi/2,   pi);
  cpfs{5} = @(x,y) cpArc(x,y, radii(3), [-1+radii(3) -1+radii(3)], pi,     3*pi/2);
  cpfs{7} = @(x,y) cpArc(x,y, radii(4), [ 1-radii(4) -1+radii(4)], 3*pi/2, 0);

  Pts = {};
  Pts{1} = [ 1-radii(1)   1];
  Pts{2} = [-1+radii(2)   1];
  Pts{3} = [-1   1-radii(2)];
  Pts{4} = [-1  -1+radii(3)];
  Pts{5} = [-1+radii(3)  -1];
  Pts{6} = [ 1-radii(4)  -1];
  Pts{7} = [ 1  -1+radii(4)];
  Pts{8} = [ 1   1-radii(1)];
  cpfs{2} = @(x,y) cpLineSegment2d(x,y, Pts{1}, Pts{2});
  cpfs{4} = @(x,y) cpLineSegment2d(x,y, Pts{3}, Pts{4});
  cpfs{6} = @(x,y) cpLineSegment2d(x,y, Pts{5}, Pts{6});
  cpfs{8} = @(x,y) cpLineSegment2d(x,y, Pts{7}, Pts{8});

  [cpx, cpy, dist] = cpCompoundObject(x,y, cpfs);

  %% mask for signed distance
  if (nargout < 3)
    sdist = dist;
  else
    if (max(radii) > 1)
      warning('signed distance calculation may not be correct for large radii');
    end

    wh1 = (( x-(1-radii(1)) ).^2  +  ( y-(1-radii(1)) ).^2 <= radii(1)^2);
    wh1 = ((x>=1-radii(1)) & (y>=1-radii(1))) .* wh1;

    wh2 = (( x-(-1+radii(2)) ).^2  +  ( y-(1-radii(2)) ).^2 <= radii(2)^2);
    wh2 = (x<=-1+radii(2) & y>=1-radii(2)) .* wh2;

    wh3 = (( x-(-1+radii(3)) ).^2  +  ( y-(-1+radii(3)) ).^2 <= radii(3)^2);
    wh3 = (x<=-1+radii(3) & y<=-1+radii(3)) .* wh3;

    wh4 = (( x-(1-radii(4)) ).^2  +  ( y-(-1+radii(4)) ).^2 <= radii(4)^2);
    wh4 = (x>=1-radii(4) & y<=-1+radii(4)) .* wh4;

    wh12 = ( x<=(1-radii(1)) & x>=(-1+radii(2)) & (y <= 1) & (y >= 0));
    wh23 = ( y<=(1-radii(2)) & y>=(-1+radii(3)) & (x <= 0) & (x >=-1));
    wh34 = ( x<=(1-radii(4)) & x>=(-1+radii(3)) & (y >=-1) & (y <= 0));
    wh41 = ( y<=(1-radii(1)) & y>=(-1+radii(4)) & (x <= 1) & (x >= 0));

    wh = wh1 | wh2 | wh3 | wh4;
    wh = wh | wh12 | wh23 | wh34 | wh41;

    sdist = (wh) .* (-dist) + (~wh) .* dist;
  end

  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
