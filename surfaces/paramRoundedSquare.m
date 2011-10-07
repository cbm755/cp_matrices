function [cpxx, cpyy, dist] = paramRoundedSquare(xx, yy, radii, cen)

  if (length(radii) == 1)
    radii = radii*[1 1 1 1];
  end
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

  keyboard

  [cpx, cpy, dist, cpCompoundObject(xx,yy, cpfs);