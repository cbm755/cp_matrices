function [cpx, cpy, dist, bdy] = cpCircleInterior(x, y, R, cen)
%CPCIRCLE  Closest Point function for a Disk.
%   [cpx, cpy, sdist] = cpDisk(x, y)
%      A unit disk centered at the origin.
%   [cpx, cpy, sdist] = cpDisk(x, y, R)
%      A disk of radius R centered at the origin.
%   [cpx, cpy, sdist] = cpDisk(x, y, R, CEN)
%      A disk of radius R centered at CEN = [xc,yc].
%


  % defaults
  if (nargin < 3)
    R = 1;
  end
  if (nargin < 4)
    cen = [0 0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);

  [th, r] = cart2pol(x, y);
  [cpx, cpy] = pol2cart(th, R);

  dist = r - R;
  bdy = dist>10*eps;
  dist(~bdy) = 0;
  cpx(~bdy) = x(~bdy);
  cpy(~bdy) = y(~bdy);
  
  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);


  % another approach:
  %r = sqrt(x.^2 + y.^2);
  %x2 = (r==0)*1.0 + (r~=0).*x;
  %r = (r==0)*1.0 + (r~=0).*r;
  %c = R ./ r;
  %cpx = c.*x2;
  %cpy = c.*y;
  %dist = sqrt( (cpx-x).^2 + (cpy-y).^2 );
