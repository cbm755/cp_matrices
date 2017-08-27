function [cpx,cpy,dist,bdy] = cpCutHoleRadius2d(x,y,cpf,hole_cen,hole_rad,cp_data, tol)
%cpCutHoleRadius2d  Cut holes in a curve by intersecting with a circle
%
%   Suppose "cpf" is closest point function or a handle to one.
%
%     [cpx, cpy, dist, bdy] = cpCutHole2d(x, y, cpf, ...
%         hole_cen, arclen_target, dx, bw, cp_data)
%
%   x, y, z: the points for which we compute the closest points.
%
%   hole_cen = [hx hy hz] gives the center of the hole.  If this is
%   not on the surface given by cpf, it will be projected onto the
%   surface (by calling cpf).
%
%   SA_target is the surface area we would like the hole to have.  The
%   algorithm tries to estimate a radius of a sphere such that the
%   intersection of the surface and the sphere has this surface srea.
%
%   Multiple holes are supported: pass each center as a row of
%   "hole_cen" and pass a vector to "SA_target".
%
%   Caveats:
%
%     * TODO: the grid in x, y must already be banded with bandwidth
%       bw*dx.  This is a restriction based on how we estimate the
%       surface area integral.  It should be fixed using Kublic and Tsai
%       or similar.
%
%     * See also cpCutHoleRadius.

% x,y coordinate matrices of any same dimension
%
% hole_cen = [hx hy] gives the center of the hole.  If this is not on
% the surface given by cpf, it will be projected onto the surface (by
% calling cpf).
%
% Multiple holes are supported: pass each center as a row of "hole_cen"
% and pass a vector to "arclen_target".
%
% See cpCutHole for caveats and additional documentation.

  if (nargin < 7)
    tol = 1e-10;
  end
  tol2 = min(1e-3, 100*tol);

  if (nargin < 6 || isempty(cp_data))
    [cpx, cpy] = cpf(x, y);
  else
    cpx = cp_data{1};
    cpy = cp_data{2};
    assert (isequal(size(x), size(cpx), size(cpy)), ...
            'cp data sizes must match coordinate inputs x, y, z')
  end
  assert (isequal(size(x), size(y)), 'coordinate inputs must be same size')

  bdy = zeros(size(cpx));
  [num_holes, temp] = size(hole_cen);
  assert (temp == 2)

  for k=1:num_holes
    xh = hole_cen(k, 1)
    yh = hole_cen(k, 2)
    % project onto surface
    [xh yh] = cpf(xh, yh)
    R_H = hole_rad(k)

    % index of all points whose cp are within the sphere
    I = ((cpx-xh).^2 + (cpy-yh).^2 < R_H^2);
    I = I(:);

    %% fixed-point iteration
    xl = cpx(I);
    yl = cpy(I);

    sdist = 1;

    while (abs(sdist) > tol | abs(dist1) > tol)
      % project onto sphere "hole"
      [xl2, yl2, sdist] = cpCircle(xl, yl, R_H, [xh yh]);

      % project onto surface
      [xll, yll, dist1] = cpf(xl2, yl2);

      % avoid projecting back and forth infinite times
      dist2=sqrt((xl-xll).^2 + (yl-yll).^2);
      dist3=sqrt((xl-xl2).^2 + (yl-yl2).^2);
      if  dist2  < tol2 & dist3 > R_H/2
        xll = xll+10^(-5)*rand; % add perturbation
        yll = yll+10^(-5)*rand;
      end

      xl=xll;
      yl=yll;
    end

    cpx(I) = xl;
    cpy(I) = yl;

    bdy = bdy | I;
  end

  dist = sqrt((cpx-x).^2 + (cpy-y).^2);
end
