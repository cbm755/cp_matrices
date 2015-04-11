function [cpx,cpy,cpz, sdist,bdy] = cpEllipsoidInterior(x,y,z, AB, cen)
 % defaults
  if (nargin < 4) || isempty(AB)
    AB = [1.25 0.75];
  end
  if (nargin < 5) || isempty(cen)
    cen = [0 0 0];
  end

  a = AB(1); b = AB(2); c = b;
  
  cpx = x; cpy = y; cpz = z; sdist = zeros(size(x));
  bdy = (x-cen(1)).^2/a^2 + (y-cen(2)).^2/b^2 + (z-cen(3)).^2/c^2 - 1 > 10*eps;
  
  [cpx(bdy), cpy(bdy), cpz(bdy), sdist(bdy)] = cpEllipsoid(x(bdy),y(bdy),z(bdy),AB,cen);
end