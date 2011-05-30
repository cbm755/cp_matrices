function [cpx,cpy,cpz, dist] = cpSphere(x,y,z, R, xc,yc,zc)
%CPSPHERE  Closest point function for a sphere.
%   [cpx,cpy,cpz, dist] = cpSphere(x,y,z, R) returns the
%   closest point and distance to (x,y,z).  If R is omitted it
%   defaults to a unit sphere, centered at the origin.
%
%   [cpx,cpy,cpz, dist] = cpSphere(x,y,z, R, xc,yc,zc) is a sphere
%   of radius R, centered at (xc,yc,zc)

  % default radius of 1
  if (nargin < 4)
    R = 1;
  end
  if (nargin == 5) | (nargin == 6)
    error('must specify all of (xc,yc,zc)');
  end
  % default center is the origin
  if (nargin < 7)
    xc = 0; yc = 0; zc = 0;
  end

  % shift to the origin
  x = x - xc;
  y = y - yc;
  z = z - zc;

  [th, phi, r] = cart2sph(x,y,z);
  [cpx,cpy,cpz] = sph2cart(th, phi, R);

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );

  % shift back
  cpx = cpx + xc;
  cpy = cpy + yc;
  cpz = cpz + zc;
