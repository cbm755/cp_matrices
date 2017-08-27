function [cpx,cpy,cpz,dist,bdy] = cpCutHole(x,y,z,cpf,hole_cen,SA_target,dx,bw,cp_data, tol)
%cpCutHole  Cut holes in a surface of prescribed surface area
%
%   Suppose "cpf" is closest point function or a handle to one.
%
%     [cpx, cpy, cpz, dist, bdy] = cpCutHole(x, y, z, cpf, ...
%         hole_cen, SA_target, dx, bw, cp_data)
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
%     * TODO: the grid in x, y, z must already be banded with bandwidth
%       bw*dx.  This is a restriction based on how we estimate the
%       surface area integral.  It should be fixed using Kublic and Tsai
%       or similar.
%
%     * See also cpCutHoleRadius.

  if (nargin < 10)
    tol = 1e-10;
  end

  if (nargin < 9 || isempty(cp_data))
    % TODO: might redesign this; for now, one has to to pass a band
    % so might as well pass the cp's as well.
    [cp_data{1:3}] = cpf(x, y, z);
  end
  cpx = cp_data{1};
  cpy = cp_data{2};
  cpz = cp_data{3};
  % TODO: if nargout is 5, get bdy from cpf instead
  %nargout(cpf)
  bdy = zeros(size(cpx));
  [num_holes, temp] = size(hole_cen);
  assert (temp == 3)

  for k=1:num_holes
    xh = hole_cen(k, 1)
    yh = hole_cen(k, 2)
    zh = hole_cen(k, 3)
    % project to ensure hole center is on the surface
    [xh yh zh] = cpf(xh, yh, zh)

    SA_wanted_k = SA_target(k);
    R_bar = sqrt(SA_wanted_k/pi);
    RR = linspace(0.5*R_bar, 1.2*R_bar, 1000);
    for i = 1:length(RR)
      I = ((cpx-xh).^2+(cpy-yh).^2+(cpz-zh).^2 < RR(i)^2);
      assert(nnz(I) > 0, 'empty sphere/surface intersection!')
      SASA(i) = sum(I(:))*(dx)^2/(2*bw);
    end
    %SASA./SA_wanted_k
    [tilde, i] = min(abs(SASA./SA_wanted_k-1));
    hole_rad(k) = RR(i);
    assert (i > 1 && i < length(RR), 'should be in the middle')
    assert (tilde < 0.1)
  end

  [cpx, cpy, cpz, dist, bdy] = cpCutHoleRadius(x, y, z, cpf, hole_cen, hole_rad, cp_data, tol);
end
