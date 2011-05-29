function [cpx, cpy, cpz, dist] = cpSphere(x, y, z, R)
%CPSPHERE  Closest point function for a sphere.

  [th, phi, r] = cart2sph(x,y,z);
  [cpx,cpy,cpz] = sph2cart(th, phi, R);

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
