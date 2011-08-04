function [x,y,z] = paramFishBowl(n, R, cen)
  warning('use paramSphereRing instead');

  % default
  if (nargin < 2)
    R = 1;
  end
  if (nargin < 3)
    cen = [0,0,0];
  end

  % multiple by two is the easiest way to get the equator right.
  [xs,ys,zs] = sphere(2*n);
  zt = 0.4;

  ind = find(zs(:,1)<zt);

  x = xs(ind,:);
  y = ys(ind,:);
  z = zs(ind,:);

  x = R*x + cen(1);
  y = R*y + cen(2);
  z = R*z + cen(3);
  