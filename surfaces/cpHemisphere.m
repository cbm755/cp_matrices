function [cpx,cpy,cpz, dist, bdy] = cpHemisphere(x,y,z, R, xc,yc,zc)
%CPHEMISPHERE  Closest point function for a hemisphere.
%   The hemisphere consists of those points with z >= 0.
%   [cpx,cpy,cpz, dist, bdy] = cpHemisphere(x,y,z)
%      A unit hemisphere with "center" (sphere center) at the origin.
%      "bdy" is non-zero for points on the boundary.
%   [cpx,cpy,cpz, dist, bdy] = cpHemisphere(x,y,z, R)
%      A radius R hemisphere centered at the origin.
%   [cpx,cpy,cpz, dist, bdy] = cpHemisphere(x,y,z, R, xc,yc,zc)
%      A radius R hemisphere centered at (xc,yc,zc).

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

  [cpx, cpy, cpz] = cpSphere(x, y, z, R);

  % should not include the bdy itself (hence not z <= 0)
  bdy = (z < 0);

  % points with z < 0 map to the z=0 plane circle of radius R
  [th, r, zp] = cart2pol(x, y, z);
  cpth = th;
  cpr = R*ones(size(th));
  cpzp = zeros(size(th));
  [cpx2, cpy2, cpz2] = pol2cart(cpth, cpr, cpzp);

  cpx(bdy) = cpx2(bdy);
  cpy(bdy) = cpy2(bdy);
  cpz(bdy) = cpz2(bdy);

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );

  % shift back to center
  cpx = cpx + xc;
  cpy = cpy + yc;
  cpz = cpz + zc;
