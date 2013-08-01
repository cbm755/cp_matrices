function [x,y,z,u,v] = paramMobiusStrip(N, Rad, Thick, cen, onlyedge)
%PARAMMOBIUSSTRIP   Parameterization of a Mobius strip
%   [x,y,z] = paramMobiusStrip(N) returns a mesh.  SURF(x,y,z) can
%   then be used to make a plot.
%
%   [x,y,z] = paramMobiusStrip(N, Radius, Width, Cen) adjusts the
%   "radius" and "width" of the mobius strip and shifts it to be
%   centered at "Cen".  Width here is measured relative to Radius.
%
%   [x,y,z] = paramMobiusStrip(N, [], [], [], true) returns only the
%   edge of the mobius strip, a codimension-2 curve that can be
%   plotted with plot3(x,y,z).

%   TODO: add an (optional) parameter to just return the edge?

  % defaults
  if (nargin < 2) || (isempty(Rad))
    Rad = 1;
  end
  if (nargin < 3) || (isempty(Thick))
    Thick = 0.35;
  end
  if (nargin < 4) || (isempty(cen))
    cen = [0 0 0];
  end
  if (nargin < 5)
    onlyedge = false;
  end

  N = max(3,N);

  if (onlyedge)
    u = 0:(4*pi/(N)):(4*pi);
    v = 1;
  else
    uu = 0:(2*pi/(N)):(2*pi);
    %darclen = 2*pi*(Rad+Thick*Rad)/N
    darclen = 2*pi*(Rad)/N;
    m = ceil(2*Rad*Thick / darclen);
    m = max(1, m);
    vv = -1:(2/m):1;
    [u,v] = meshgrid(uu,vv);
  end

  x = Rad*(1 + Thick*v .* cos(.5*u)) .* cos(u) - Rad*Thick/2;
  y = Rad*(1 + Thick*v .* cos(.5*u)) .* sin(u);
  z = Rad*Thick*v .* sin(0.5*u);

  x = x + cen(1);
  y = y + cen(2);
  z = z + cen(3);

