function [cpx,cpy,cpz, dist] = cpSimpleSphere(x,y,z, R)
%CPSIMPLESPHERE  Closest point function for a sphere.
%   [cpx,cpy,cpz, dist] = cpSphere(x,y,z, R) returns the
%   closest point and distance to (x,y,z).  If R is omitted it
%   defaults to a unit sphere, centered at the origin.

  % default radius of 1
  if (nargin < 4)
    R = 1;
  end

  [th, phi, r] = cart2sph(x,y,z);
  [cpx,cpy,cpz] = sph2cart(th, phi, R);

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
