function [cpx,cpy,sdist,bdy] = cpEllipseInterior(x,y, AB, cen)
 % defaults
  if (nargin < 3) || isempty(AB)
    AB = [1.25 0.75];
  end
  if (nargin < 4) || isempty(cen)
    cen = [0 0];
  end

  a = AB(1); b = AB(2);
  
  cpx = x; cpy = y;  sdist = zeros(size(x));
  bdy = (x-cen(1)).^2/a^2 + (y-cen(2)).^2/b^2 - 1 > 10*eps;
  
  [cpx(bdy), cpy(bdy), sdist(bdy)] = cpEllipse(x(bdy),y(bdy),a,b,cen);
end