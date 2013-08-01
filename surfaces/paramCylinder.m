function [x,y,z] = paramCylinder(n, zlim, R, cen)
%PARAMCYLINDER   A parameterization of a cylinder

  % default radius
  if (nargin < 3),  R = 1;  end
  % default bottom/top
  if (nargin < 2),  zlim = [-1  1];  end
  % default center (in x,y) is the origin
  if (nargin < 4),  cen = [0,0];  end

  zlo = zlim(1);  zhi = zlim(2);

  theta = linspace (0, 2*pi, n+1);
  % estimate how many grids to use in the vertical direction
  sidelen = 2*pi*R/n;
  nv = ceil((zhi - zlo) / sidelen);
  z1d = linspace (zlo, zhi, nv+1);

  [theta, z] = meshgrid(theta, z1d);

  x = R .* cos(theta);
  y = R .* sin(theta);

  x = x + cen(1);
  y = y + cen(2);
