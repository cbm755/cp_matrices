function [x,y,z] = paramEllipsoid(N, AB, cen, ax)
%PARAMELLIPSOID  A parameterization of an ellipsoid
%   [x,y,z] = paramEllipsoid(N, [a b]) returns a mesh for a ellipsoid
%   major and minor principal axes.  If omitted [a b] take default
%   values.  surf(x,y,z) can be used to make a plot.
%
%   [x,y,z] = paramEllipsoid(N, [a b], CEN) returns a mesh centered at CEN.
%

  % defaults
  if (nargin < 2) || isempty(AB)
    AB = [1.25 0.75];
  end
  if (nargin < 3) || isempty(cen)
    cen = [0 0 0];
  end
  if (nargin < 4)
    ax = 'x';
  end

  a = AB(1);
  b = AB(2);

  % TODO: not visually pleasing when the ratio a/b is small
  if (a > b)
    M = ceil(b/a*N*2);
  else
    M = N;
    N = ceil(a/b*M/2);
  end
  N = max(3,N);
  M = max(3,M);

  % TODO: some kind of tighter-than-chebychev spacing, might look
  % better---logspace in th?
  th = 0:(2*pi/M):(2*pi);
  u = 0:(pi/N):pi;
  [u2,th2] = meshgrid(u,th);

  switch ax
    case 'x'
      x = a*cos(u2);
      y = b*sin(u2) .* cos(th2);
      z = b*sin(u2) .* sin(th2);
    case 'y'
      y = a*cos(u2);
      x = b*sin(u2) .* cos(th2);
      z = b*sin(u2) .* sin(th2);
    case 'z'
      z = a*cos(u2);
      y = b*sin(u2) .* cos(th2);
      x = b*sin(u2) .* sin(th2);
    otherwise
      error('invalid axis of revolution');
  end

  x = x + cen(1);
  y = y + cen(2);
  z = z + cen(3);
