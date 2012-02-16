function [cpx,cpy,cpz, sdist] = cpSphere(x,y,z, R, cen)
%CPSPHERE  Closest point function for a sphere.
%   [cpx,cpy,cpz, sdist] = cpSphere(x,y,z, R) returns the
%   closest point and distance to (x,y,z).  If R is omitted it
%   defaults to a unit sphere, centered at the origin.
%
%   [cpx,cpy,cpz, sdist] = cpSphere(x,y,z, R, [xc,yc,zc]) is a sphere
%   of radius R, centered at (xc,yc,zc)
%
%   Note: returns signed distance (with negative inside).


  % defaults
  if (nargin < 4)
    R = 1;
  end
  if (nargin < 5)
    cen = [0, 0, 0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  [th, phi, r] = cart2sph(x,y,z);
  [cpx,cpy,cpz] = sph2cart(th, phi, R);

  %dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
  sdist = r - R;

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);
