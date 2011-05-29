function [cpx, cpy, dist] = cpCircle(x, y, R)
%CPCIRCLE  Closest Point function for a circle.

  [th, r] = cart2pol(x, y);
  [cpx, cpy] = pol2cart(th, R);

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 );

  % another approach:
  %r = sqrt(x.^2 + y.^2);
  %x2 = (r==0)*1.0 + (r~=0).*x;
  %r = (r==0)*1.0 + (r~=0).*r;
  %c = R ./ r;
  %cpx = c.*x2;
  %cpy = c.*y;
  %dist = sqrt( (cpx-x).^2 + (cpy-y).^2 );

