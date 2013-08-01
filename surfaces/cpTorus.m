function [cpx,cpy,cpz, sdist] = cpTorus(x,y,z, A, B, cen)
%CPTORUS  Closest point function for a torus
%   [cpx,cpy,cpz, sdist] = cpTorus(x,y,z, R, r) returns the closest
%   point and distance to (x,y,z) for a torus with major radius R and
%   minor radius r.  R and r default to 1 and 0.4 respectively if
%   omitted.
%
%   [cpx,cpy,cpz, sdist] = cpTorus(x,y,z, R, r, [xc yc zc]) is as
%   above but centered at (xc,yc,zc).
%
%   Note: returns signed distance (with negative inside).

  % defaults
  if (nargin < 4)
    A = 1;
  end
  if (nargin < 5)
    B = 0.4;
  end
  if (nargin < 6)
    cen = [0 0 0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  % cyclindrical coors
  [th, r, tilde] = cart2pol(x,y,z);

  % now we can work in a 2D plane
  [cpr, cpz, sdist] = cpCircle(r, z, B, [A 0]);

  % and convert back
  [cpx,cpy] = pol2cart(th, cpr);

  %dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
  %assert(all(all(all(dist - abs(sdist) <= 1e-15))))

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);
