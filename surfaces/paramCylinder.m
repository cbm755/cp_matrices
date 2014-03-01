function [x,y,z] = paramCylinder(n, zlim, R, cen)
%PARAMCYLINDER   A parameterization of a cylinder

  % default radius
  if (nargin < 3),  R = 1;  end
  % default bottom/top
  if (nargin < 2),  zlim = [-1  1];  end
  % default center (in x,y) is the origin
  if (nargin < 4),  cen = [0,0];  end

  zlo = zlim(1);  zhi = zlim(2);

  % depending on which way is larger
  if (2*pi*R > (zhi-zlo))
    nt = n;
    % estimate how many grids in the vertical direction
    nz = ceil((zhi-zlo) / (2*pi*R/nt));
  else
    nz = n;
    % estimate how many grids in the angular direction
    nt = ceil(2*pi*R / ((zhi-zlo)/nz));
  end
  nz = max(nz,4);
  z1d = linspace (zlo, zhi, nz+1);
  nt = max(nt,8);
  theta = linspace (0, 2*pi, nt+1);

  [theta, z] = meshgrid(theta, z1d);

  x = R .* cos(theta);
  y = R .* sin(theta);

  x = x + cen(1);
  y = y + cen(2);
