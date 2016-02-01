function [cpx,cpy,cpz,dist,bdy] = cpSphereRing(x,y,z,zlim,R,cen)
%CPSPHERERING   Closest point fcn for a ring, subset of sphere
%   This can be used to make Rings, Fishbowls and other things
%   that are subsets of a sphere
%
%   Code is vectorized: any size/shape for x should work.


  % defaults
  if (nargin < 5)
    R = 1;
  end
  if (nargin < 6)
    cen = [0,0,0];
  end
  if (nargin < 4)
    % default is a fishbowl
    zlim = [-2*R  0.5];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);
  zlim = zlim - cen(3);

  [cpx, cpy, cpz] = cpSphere(x, y, z, R);

  %zlo = R*sin(phis(1));
  %zhi = R*sin(phis(2));
  %radlo = R*cos(phis(1));
  %radhi = R*cos(phis(2));

  % should not include the bdy itself (hence <=)
  bdy1 = (cpz < zlim(1));
  bdy2 = (cpz > zlim(2));

  philo = asin(max(-1, zlim(1)/R));
  phihi = asin(min(1, zlim(2)/R));
  philim = [philo  phihi];
  radlim = [cos(philim(1))  cos(philim(2))];

  [th, tilde, tilde] = cart2pol(x, y, z);
  cpth = th;
  cpr = radlim(1)*ones(size(th));
  cpzp = zlim(1)*ones(size(th));
  [cpx1, cpy1, cpz1] = pol2cart(cpth, cpr, cpzp);

  [th, tilde, tilde] = cart2pol(x, y, z);
  cpth = th;
  cpr = radlim(2)*ones(size(th));
  cpzp = zlim(2)*ones(size(th));
  [cpx2, cpy2, cpz2] = pol2cart(cpth, cpr, cpzp);

  cpx(bdy1) = cpx1(bdy1);
  cpy(bdy1) = cpy1(bdy1);
  cpz(bdy1) = cpz1(bdy1);

  cpx(bdy2) = cpx2(bdy2);
  cpy(bdy2) = cpy2(bdy2);
  cpz(bdy2) = cpz2(bdy2);

  bdy = (bdy1 | bdy2);

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );

  % shift back to center
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);
