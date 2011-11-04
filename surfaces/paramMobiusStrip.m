function [x,y,z,varargout] = paramMobiusStrip(N, Rad, Thick, cen)
%PARAMMOBIUSSTRIP   Parameterization of a Mobius strip
%   [x,y,z] = paramMobiusStrip(N) returns a mesh.  SURF(x,y,z) can
%   then be used to make a plot.
%
%   [x,y,z] = paramMobiusStrip(N, Radius, Width) adjusts the "radius" and
%   "width' of the mobius strip.  Width here is measured relative to
%   Radius.

  % defaults
  if (nargin < 2)
    Rad = 1;
  end
  if (nargin < 3)
    Thick = 0.35;
  end
  if (nargin < 4)
    cen = [0 0 0];
  end

  N = max(3,N);

  uu = 0:(2*pi/(N)):(2*pi);
  %darclen = 2*pi*(Rad+Thick*Rad)/N
  darclen = 2*pi*(Rad)/N;
  m = ceil(2*Rad*Thick / darclen);
  m = max(1, m);

  vv = -1:(2/m):1;


  [u,v] = meshgrid(uu,vv);

  x = Rad*(1 + Thick*v .* cos(.5*u)) .* cos(u) - Rad*Thick/2;
  y = Rad*(1 + Thick*v .* cos(.5*u)) .* sin(u);
  z = Rad*Thick*v .* sin(0.5*u);

  x = x + cen(1);
  y = y + cen(2);
  z = z + cen(3);

  % optionally output the parameters too
  if nargout == 5
    varargout{1} = u;
    varargout{2} = v;
  end