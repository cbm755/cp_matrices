function [x,y,z] = paramVase(N, lim, ab, cen)
%PARAMVASE  A parameterization of a vase
%   [x,y,z] = paramVase(N, [th1 th2], [a b]) returns a mesh for a
%   vase, see cpVase for meaning of parameters.
%
%   surf(x,y,z) can be used to make a plot.

  % defaults
  if (nargin < 2) || isempty(lim)
    lim = [-.7*pi/2 .4*pi/2];
  end
  if (nargin < 3) || isempty(ab)
    ab = [1.5 0.75];
  end
  if (nargin < 4) || isempty(cen)
    cen = [0 0 -ab(1)/2*sum(sin(lim))];
  end

  a = ab(1);
  b = ab(2);

  % TODO: not visually pleasing when the ratio a/b is small
  if (a > b)
    M = ceil(b/a*N*2);
  else
    M = N;
    N = ceil(a/b*M/2);
  end
  N = max(3,N);
  M = max(3,M);

  th = 0:(2*pi/M):(2*pi);
  u = lim(1):((lim(2)-lim(1))/M):lim(2);
  [u2,th2] = meshgrid(u,th);

  z = a*sin(u2);
  x = b*cos(u2) .* cos(th2);
  y = b*cos(u2) .* sin(th2);

  x = x + cen(1);
  y = y + cen(2);
  z = z + cen(3);

