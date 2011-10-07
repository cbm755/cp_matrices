function [cpxx, cpyy, dist] = cpRoundedSquare(x, y, radii, cen)
%CPROUNDEDSQUARE  Closest Point function for a square with rounded corners
%   [cpx, cpy, dist] = cpRoundedSquare(x, y, radius)
%      A square with side length 2 centered centered at
%      the origin, each corner is rounded with 'radius'
%   [cpx, cpy, dist] = cpSquare(x, y, radii)
%      A square with side length 2 centered centered at
%      the origin, the corner of the jth quadrant is rounded with
%      radius 'radii(j)'.
%   [cpx, cpy, dist] = cpSquare(x, y, radii, [xc yc])
%      A square with side length 2 centered centered at
%      the point (xc,yc)

  if (length(radii) == 1)
    radii = radii*[1 1 1 1];
  end

  x = x - cen(1);
  y = y - cen(2);

  cpfs = {};
  cpfs{1} = @(x,y) cpArc(x,y, radii(1), [ 1  1] - radii(1), 0, pi/2);
  cpfs{3} = @(x,y) cpArc(x,y, radii(2), [-1  1] - radii(2), pi/2, pi);
  cpfs{5} = @(x,y) cpArc(x,y, radii(3), [-1 -1] - radii(3), pi, 3*pi/2);
  cpfs{7} = @(x,y) cpArc(x,y, radii(4), [ 1 -1] - radii(4), 3*pi/2, 0);

  Pts = {};
  Pts{1} = [ 1-radii(1)   1];
  Pts{2} = [-1+radii(2)   1];
  Pts{3} = [-1   1-radii(2)];
  Pts{4} = [-1  -1+radii(3)];
  Pts{1} = [-1+radii(3)  -1];
  Pts{2} = [ 1-radii(4)  -1];
  Pts{3} = [ 1  -1+radii(4)];
  Pts{4} = [ 1   1-radii(1)];
  cpfs{2} = @(x,y) cpLineSegment(x,y, Pts{1}, Pts{2});
  cpfs{4} = @(x,y) cpLineSegment(x,y, Pts{3}, Pts{4});
  cpfs{6} = @(x,y) cpLineSegment(x,y, Pts{5}, Pts{6});
  cpfs{8} = @(x,y) cpLineSegment(x,y, Pts{7}, Pts{8});

  [cpx, cpy, dist] = CompoundObject(x,y, cpfs);

  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
