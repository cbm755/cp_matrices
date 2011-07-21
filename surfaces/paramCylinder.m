function [x,y,z] = paramCylinder(n, R, zlo, zhi,cen)
%PARAMCYLINDER   A parameterization of a cylinder

  % default radius
  if (nargin < 2)   R = 1;   end
  % default bottom
  if (nargin < 3)   zlo = -1;   end
  % default top
  if (nargin < 4)   zhi = 1;   end
  % default center
  if (nargin < 5)   cen = [0,0];   end

  theta = linspace (0, 2*pi, n+1);
  % estimate how many grids to use in the vertical direction
  sidelen = 2*pi*R/n;
  nv = ceil((zhi - zlo) / sidelen);
  z1d = linspace (zlo, zhi, nv+1);

  [theta, z] = meshgrid (theta, z1d);

  x = R .* cos (theta);
  y = R .* sin (theta);
  %z = z;

  x = x + cen(1);
  y = y + cen(2);
