function [cpx,cpy,cpz, dist, bdy] = cpHemisphere(x,y,z, R)
%CPHEMISPHERE  Closest point function for a hemisphere.
%   The semicircle consists of those points with y >= 0.
%   "bdy" is non-zero for points on the boundary.
%   TODO: could use varargin and default to R=1 and also handle
%   other centers.

% shift to origin (if centered elsewhere)
%x = x - 2;
%y = y - 1;
%z = z - 0.5;

[cpx, cpy, cpz] = cpSphere(x, y, z, R);

% should not include the bdy itself (hence not z <= 0)
bdy = (z < 0);

[th, r, zp] = cart2pol(x, y, z);
cpth = th;
cpr = R*ones(size(th));
cpzp = zeros(size(th));
[cpx2, cpy2, cpz2] = pol2cart(cpth, cpr, cpzp);

cpx(bdy) = cpx2(bdy);
cpy(bdy) = cpy2(bdy);
cpz(bdy) = cpz2(bdy);

dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );

% shift back to center (if not origin)
%cpx = cpx + 2;
%cpy = cpy + 1;
%cpz = cpz + 0.5;
