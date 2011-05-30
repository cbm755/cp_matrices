function [cpx,cpy,cpz, dist] = cpSphere(x,y,z, varargin)
%CPSPHERE  Closest point function for a sphere.
%   [cpx,cpy,cpz, dist] = cpSphere(x,y,z, R) returns the
%   closest point and distance to (x,y,z).  If R is omitted it
%   defaults to a unit sphere, centered at the origin.
%
%   [cpx,cpy,cpz, dist] = cpSphere(x,y,z, R, xc,yc,zc) is a sphere
%   of radius R, centered at (xc,yc,zc)

  % value for R, or default to 1
  if (nargin >= 4)
    R = varargin{1};
  else
    R = 1;
  end

  % center shift, todo
  if (nargin >= 5)
    error('todo');
    xc = varargin{2}
    yc = varargin{3}
    zc = varargin{4}
  else
    xc = 0; yc = 0; zc = 0;
  end

  %x = x - xc;
  %y = y - yc;
  %z = z - zc;

  [th, phi, r] = cart2sph(x,y,z);
  [cpx,cpy,cpz] = sph2cart(th, phi, R);

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
