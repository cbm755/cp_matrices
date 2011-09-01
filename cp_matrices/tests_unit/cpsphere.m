function [cpx, cpy, cpz, dist] = cpsphere(x, y, z, R)
%
  [th,phi,r] = cart2sph(x,y,z);
  [cpx,cpy,cpz] = sph2cart(th,phi,R);

  dist = sqrt( (cpx-x).^2 + (cpy-y).^2 + (cpz-z).^2 );

