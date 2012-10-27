function [x,y,z] = paramTorus(N, R, r, cen)
%PARAMTORUS  A parameterization of a torus
%   [x,y,z] = paramTorus(N, R, r) returns a mesh for a torus with
%   major radius R and minor radius r.  R and r default to 1 and 0.4
%   respectively if omitted.  surf(x,y,z) can be used to make a plot.
%
%   [x,y,z] = paramTorus(N, R, r, cen) returns a mesh centered at cen.

  if (nargin < 2)
    R = 1;
  end
  if (nargin < 3)
    r = 0.4;
  end
  if (nargin < 4)
    cen = [0, 0, 0];
  end

  M = N;              % points around major radius
  N = ceil(N/(R/r));  % points around minor radius
  M = max(8, M);      % minimums of each
  N = max(8, N);

  x = zeros(N+1,M+1);
  y = x;
  z = x;

  dth = 2*pi/N;
  for i=0:N
    th = i*dth;

    x(i+1,:) = (R+r*sin(th))*cos((0:M)*2*pi/M);
    y(i+1,:) = (R+r*sin(th))*sin((0:M)*2*pi/M);
    z(i+1,:) = -r*cos(th)*ones(1,M+1);
  end

  x = x + cen(1);
  y = y + cen(2);
  z = z + cen(3);

