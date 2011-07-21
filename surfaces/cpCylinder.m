function [cpx,cpy,cpz, dist, bdy] = cpCylinder(x,y,z, R, zlo, zhi, cen)
%CPCYLINDER  Closest point function for a cylinder without caps
%   ...
%   TODO: version with end-caps?
%
%   Code is vectorized: any size/shape for x should work.


  % default radius
  if (nargin < 4)   R = 1;   end
  % default bottom
  if (nargin < 5)   zlo = -1;   end
  % default top
  if (nargin < 6)   zhi = 1;   end
  % default center (in x,y) is the origin
  if (nargin < 7)   cen = [0,0];   end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  %z = z;

  bdy1 = (z < zlo);
  bdy2 = (z > zhi);

  [th, r, zp] = cart2pol(x, y, z);
  cpth = th;
  cpr = R*ones(size(th));
  cpzp = zp;
  cpzp(bdy1) = zlo;
  cpzp(bdy2) = zhi;

  [cpx, cpy, cpz] = pol2cart(cpth, cpr, cpzp);

  bdy = bdy1 | bdy2;

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );

  % shift back to center
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  %cpz = cpz;
