function [cpx, cpy, dist] = cpPoint2d(x, y, cen)
%CPPOINT2D  Closest Point function for a point
%
%   Code is vectorized: any size/shape for x should work.


  % defaults
  if (nargin < 3)
    cen = [0,0];
  end

  cpx = cen(1)*ones(size(x));
  cpy = cen(2)*ones(size(x));

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 );
