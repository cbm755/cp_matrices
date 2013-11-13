function [cpx,cpy,cpz,dist,bdy] = cpVase(x, y, z, lim, ab, cen)
%CPVASE   Closest point function for a vase
%   The vase is a subset of an oblate/prolate ellipsoid.
%   [cpx,cpy,cpz,dist,bdy] = cpVase(x, y, z, [th1 th2], [a b])
%      part of a ellipsoid with axis 'a' along 'z' and axis 'b' in
%      the x-y plane.  The ellipse is truncated by [th1 th2] in
%      the range [-pi/2, pi/2], specifically:
%          bottom of vase is a*sin(th1)
%          top of vase is a*sin(th2)
%   [cpx,cpy,cpz,dist,bdy] = cpVase(..., cen)
%      same but centered at 'cen', center is measured from the center
%      of the ellipsoid.  If omitted, by default it centers the vase.
%
%   Code is vectorized: any size/shape for x should work.

  % defaults
  if (nargin < 4) || isempty(lim)
    lim = [-.7*pi/2 .4*pi/2];
  end
  if (nargin < 5) || isempty(ab)
    ab = [1.5 0.75];
  end
  if (nargin < 6) || isempty(cen)
    cen = [0 0 -ab(1)/2*sum(sin(lim))];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  %this is basically cpSurfOfRevolution()
  cpf = @cpEllipseCurve;
  [th,r,zz] = cart2pol(x,y,z);
  % sin/cos switch here (cpEllipseCurve uses x axis)
  lim = pi/2 - [lim(2) lim(1)];
  [cpzz, cpr, dist, bdy] = cpf(zz, r, lim, ab);
  [cpx,cpy,cpz] = pol2cart(th, cpr, cpzz);

  % shift back to center
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);

