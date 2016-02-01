function [x,y,z] = paramSphereRing(n, zlim, R, cen)
%PARAMSPHEREIRING   Parameterization of part of a sphere
%   [x,y,z] = paramSphereRing(N) returns a mesh for fishbowl of unit
%   radius.  SURF(x,y,z) can be used to make a plot.
%
%   [x,y,z] = paramSphereRing(N, ZLIM, R, CEN) returns part of a
%   spherical mesh where sphere has radius R and is centered at the
%   3-vector CEN.  The section of the sphere extends from ZLIM(1) to
%   ZLIM(2).  ZLIM defaults to [-2*R 0.5], the fishbowl.  Note its
%   safer (?) to give one limit as -2*R rather than -R.

  % defaults
  if (nargin < 3)
    R = 1;
  end
  if (nargin < 4)
    cen = [0 0 0];
  end
  if (nargin < 2)
    % default is a fishbowl
    zlim = [-2*R  0.5];
  end

  % shift to origin as zlim is in the actual coordinates
  zlim = zlim - cen(3);

  philo = asin(max(-1,zlim(1)/R));
  phihi = asin(min(1,zlim(2)/R));
  philim = [philo phihi];

  th = linspace(0, 2*pi, n+1);

  % how many divisions in the phi (vertical) direction to match the
  % n divisions in the horizontal
  dth = th(2) - th(1);
  m = ceil((philim(2) - philim(1)) / dth);
  phi = linspace(philim(1), philim(2), m+1);

  x = R * (cos(th)' * cos(phi)) + cen(1);
  y = R * (sin(th)' * cos(phi)) + cen(2);
  z = R * (ones(size(th))' * sin(phi)) + cen(3);
