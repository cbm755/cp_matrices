function [cpx,cpy,cpz, dist, bdy] = cpCylinder(x,y,z, zlim, R, cen)
%CPCYLINDER  Closest point function for a cylinder without caps
%   A cylinder rising in the z-direction.
%   [cpx,cpy,cpz, dist, bdy] = cpCylinder(x,y,z, zlim, R, cen)
%     'zlim' defaults to [-1 1]
%     radius 'R' defaults to 1.
%     'cen', location in xy-plane of cylinder, default: [0,0].
%
%   TODO: version with end-caps?
%
%   Code is vectorized: any size/shape for x should work.


  % default radius
  if (nargin < 5),   R = 1;   end
  % default bottom/top
  if (nargin < 4),   zlim = [-1  1];   end
  % default center (in x,y) is the origin
  if (nargin < 6),   cen = [0,0];   end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);

  zlo = zlim(1);
  zhi = zlim(2);

  bdy1 = (z < zlo);
  bdy2 = (z > zhi);

  [th, tilde, zp] = cart2pol(x, y, z);
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
