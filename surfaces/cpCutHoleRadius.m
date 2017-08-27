function [cpx,cpy,cpz,dist,bdy] = cpCutHoleRadius(x,y,z,cpf,hole_cen,hole_rad,cp_data,tol)
%cpCutHole  Cut holes in a surface by intersecting with a sphere
%
%   See "cpCutHole" if you'd prefer to specify surface area of the hole.
%
%   Suppose "cpf" is closest point function or a handle to one.
%
%     [cpx, cpy, cpz, dist, bdy] = cpCutHoleRadius(x, y, z, cpf, ...
%           hole_cen, hole_rad)
%
%     [...] = cpCutHoleRadius(..., hole_rad, cp_data)
%
%     [...] = cpCutHoleRadius(..., hole_rad, cp_data, tol)
%
%   x, y, z: the points for which we compute the closest points.
%
%   hole_cen = [hx hy hz] gives the center of the hole.  If this is
%   not on the surface given by cpf, it will be projected onto the
%   surface (by calling cpf).
%
%   hole_rad is the radius of a sphere we place at "hole_cen".  The
%   intersection of the surface and this sphere will be cut.
%
%   Multiple holes are supported: pass each center as a row of
%   "hole_cen" and pass a vector to "hole_rad".  However, they should not
%   self-intersect.  TODO: should be easy enough to check this.
%
%   Existing arrays cpx, cpy and cpz can be passed as a cell array
%   "cp_data = {cpx, cpy, cpz}".  If "cp_data" is omitted, this data will
%   be calculated by calling "cpf".
%
%   Caveats:
%
%     * Our algorithm is based on a fixed point iteration between a
%       sphere and the surface.  The two must intersect reasonably close
%       to orthogonal for convergence.
%
%     * TODO: we lose the original boundary information if the original
%       surface was open.
%
%     * TODO: multiple holes should be labeled with different integers
%       in "bdy".

  if (nargin < 8)
    tol = 1e-10;
  end
  tol2 = min(1e-3, 100*tol);

  if (nargin < 7 || isempty(cp_data))
    [cpx, cpy, cpz] = cpf(x, y, z);
  else
    cpx = cp_data{1};
    cpy = cp_data{2};
    cpz = cp_data{3};
    assert (isequal(size(x), size(cpx), size(cpy), size(cpz)), ...
            'cp data sizes must match coordinate inputs x, y, z')
  end
  assert (isequal(size(x), size(y), size(z)), 'coordinate inputs must be same size')

  bdy = zeros(size(cpx));
  [num_holes, temp] = size(hole_cen);
  assert (temp == 3)

  for k=1:num_holes
    xh = hole_cen(k, 1);
    yh = hole_cen(k, 2);
    zh = hole_cen(k, 3);
    R_H = hole_rad(k);
    % project to ensure hole center is on the surface
    [xh yh zh] = cpf(xh, yh, zh)

    % index of all points whose cp are within the sphere
    I = ((cpx-xh).^2 + (cpy-yh).^2 + (cpz-zh).^2 < R_H^2);
    I = I(:);

    % l for loop, points to project onto boundary of hole on surface
    xl = cpx(I);
    yl = cpy(I);
    zl = cpz(I);
    
    
    sdist = 1;
    
    while  abs(sdist) > tol | abs(dist1) > tol
       
        % project onto sphere "hole"
        [xl2,yl2,zl2,sdist] = cpSphere(xl,yl,zl,R_H,[xh,yh,zh]);
                        
        % project onto surface
        [xll,yll,zll,dist1] = cpf(xl2, yl2, zl2);

        % avoid projecting back and forth infinite times
        dist2=sqrt((xl-xll).^2 + (yl-yll).^2 + (zl-zll).^2);
        dist3=sqrt((xl-xl2).^2 + (yl-yl2).^2 + (zl-zl2).^2);
        if  dist2  < tol2 & dist3 > R_H/2
           xll = xll+10^(-5)*rand; % add perturbation
           yll = yll+10^(-5)*rand;
           zll = zll+10^(-5)*rand;

        end
        
        xl=xll;
        yl=yll;
        zl=zll; 
          
    end
       
    cpx(I) = xl;
    cpy(I) = yl;
    cpz(I) = zl;
    
    bdy = bdy | I; 
  end

  dist = sqrt((cpx-x).^2 + (cpy-y).^2 + (cpz-z).^2);
end
