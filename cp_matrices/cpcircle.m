function [cpx, cpy, dist] = cpcircle(x, y, R)
%
  [th,r] = cart2pol(x,y);
  [cpx,cpy] = pol2cart(th, r);

  dist = sqrt( (cpx-x).^2 + (cpy-y).^2 );

