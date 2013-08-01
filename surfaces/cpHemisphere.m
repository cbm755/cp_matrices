function [cpx,cpy,cpz, dist, bdy] = cpHemisphere(x,y,z, R, cen)
%CPHEMISPHERE  Closest point function for a hemisphere.
%   The hemisphere consists of those points with z >= 0.
%   [cpx,cpy,cpz, dist, bdy] = cpHemisphere(x,y,z)
%      A unit hemisphere with "center" (sphere center) at the origin.
%      "bdy" is non-zero for points on the boundary.
%   [cpx,cpy,cpz, dist, bdy] = cpHemisphere(x,y,z, R)
%      A radius R hemisphere centered at the origin.
%   [cpx,cpy,cpz, dist, bdy] = cpHemisphere(x,y,z, R, CEN)
%      A radius R hemisphere centered at CEN.
%
%   Code is vectorized: any size/shape for x should work.


  % default radius of 1
  if (nargin < 4)
    R = 1;
  end
  % default center is the origin
  if (nargin < 5)
    cen = [0,0,0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);


  [cpx, cpy, cpz] = cpSphere(x, y, z, R);

  % should not include the bdy itself (hence <=)
  bdy = (z < 0);

  % points with z < 0 map to the z=0 plane circle of radius R
  [th, tilde, tilde] = cart2pol(x, y, z);
  cpth = th;
  cpr = R*ones(size(th));
  cpzp = zeros(size(th));
  [cpx2, cpy2, cpz2] = pol2cart(cpth, cpr, cpzp);

  cpx(bdy) = cpx2(bdy);
  cpy(bdy) = cpy2(bdy);
  cpz(bdy) = cpz2(bdy);

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );

  % shift back to center
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);
