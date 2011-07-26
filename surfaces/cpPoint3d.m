function [cpx,cpy,cpz, dist] = cpPoint3d(x,y,z, cen)
%CPPOINT3D  Closest Point function for a single point
%   Code is vectorized: any size/shape for x should work.


  % defaults
  if (nargin < 4)
    cen = [0,0,0];
  end

  cpx = cen(1)*ones(size(x));
  cpy = cen(2)*ones(size(x));
  cpz = cen(3)*ones(size(x));

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
