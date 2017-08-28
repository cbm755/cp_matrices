function [cpx,cpy,cpz,dist,bdy] = cpCutPlane(xx, yy, zz, cpf, point, normal, cp_data)
%cpCutPlane  Cut a surface by a oriented plane
%
%   Suppose "cpf" is closest point function or a handle to one.
%
%     [cpx, cpy, cpz, dist, bdy] = cpCutPlane(x, y, z, cpf)
%
%     [...] = cpCutPlane(x, y, z, cpf, plane_pt)
%
%     [...] = cpCutPlane(x, y, z, cpf, plane_pt, plane_normal)
%
%     [...] = cpCutPlane(x, y, z, cpf, plane_pt, plane_normal, cpdata)
%
%   x, y, z: the points for which we compute the closest points.
%
%   plane_pt: a point that the plane goes through.  Defaults to the
%   origin if omitted or empty.
%
%   plane_normal: normal vector to the plane.  Defaults to "[0; 0; 1]"
%   if empty or omitted.  All parts of the surface below the plane
%   will be cut ("below" means in the opposite direction as the normal).
%
%   cp_data:  Existing arrays cpx, cpy and cpz can be passed as a cell
%   array "cp_data = {cpx, cpy, cpz}".  If "cp_data" is omitted, it will
%   be calculated by calling "cpf".
%
%   Caveats:
%
%     * Our algorithm is based on a fixed point iteration between the
%       surface and the plane.  The two must intersect reasonably close
%       to orthogonal for convergence.
%
%     * TODO: we lose the original boundary information if the original
%       surface was open.
%
%     * TODO: normal harcoded to [0 0 1].

  tol = 1e-14;
  tol2 = min(1e-3, 1000*tol);

  % draw the points as they converge (uses figure(1))
  make_plots = true;

  if (nargin < 5 || isempty(point))
    point = [0; 0; 0];
  end

  if (nargin < 6 || isempty(normal))
    normal = [0; 0; 1];
  else
    warning('WIP: normal hardcoded to 0 0 1')
  end

  if (nargin < 7 || isempty(cp_data))
    [cpx, cpy, cpz] = cpf(xx, yy, zz);
  else
    cpx = cp_data{1};
    cpy = cp_data{2};
    cpz = cp_data{3};
    assert (isequal(size(xx), size(cpx), size(cpy), size(cpz)), ...
            'cp data sizes must match coordinate inputs x, y, z')
  end
  assert (isequal(size(xx), size(yy), size(zz)), 'coordinate inputs must be same size')

  bdy = zeros(size(cpx));  % TODO: assumed zero, should accept in cp_data

  assert (numel(point) == 3)
  assert (numel(normal) == 3)

  xh = point(1);  yh = point(2);  zh = point(3);

  % index of all points whose cp are below the plane
  bdy = (cpz - zh) < 0;
  bdy = bdy(:);

  cx = cpx(bdy);
  cy = cpy(bdy);
  cz = cpz(bdy);

  totalpts = nnz(bdy);

  if (make_plots)
    figure(1);
    hh1 = plot3(cx,cy,cz, 'rx');
    hh2 = plot3(cx,cy,cz, 'k.');
  end


  iter = 0;
  all_dists = inf*ones(size(cx));
  all_distp = all_dists;

  while (true)
    % points that have not yet converged
    I = abs(all_dists) > tol | abs(all_distp) > tol;
    % points we're curently working on, currently on the surface
    x = cx(I);
    y = cy(I);
    z = cz(I);

    if (make_plots)
      set(hh1, 'xdata', x, 'ydata', y, 'zdata', z)
      set(hh2, 'xdata', cx(~I), 'ydata', cy(~I), 'zdata', cz(~I))
      drawnow
    end

    disp(sprintf('iter%2d: %d of %d points have converged, norms = [%.2g, %.2g]', ...
         iter, nnz(~I), totalpts, norm(all_dists(I), 2), norm(all_distp(I), 2)));
    if (nnz(I) == 0)
      break
    end
    iter = iter + 1;

    % project onto plane
    %[x1, y1, z1, sdist] = cpPlane(x, y, z, point, normal);
    x1 = x;
    y1 = y;
    distp = abs(z - zh);
    z1 = ones(size(z))*zh;

    % project onto surface
    [x2, y2, z2, dists] = cpf(x1, y1, z1);

    % cp for surface and cp for plane should be converging to each other
    % To avoid (rare?) projecting back and forth infinitly, we check and
    % perturb.
    dist2=sqrt((x-x2).^2 + (y-y2).^2 + (z-z2).^2);
    if (any(dist2 < tol2 & distp > 0.1))  % TODO: should depend on dx?
      warning('might not be converged, do we handle this well enough?')
      II = (dist2 < tol2 & distp > 0.1);
      x2(II) = x2(II) + 1e-5*rand(nnz(II),1); % add perturbation
      y2(II) = y2(II) + 1e-5*rand(nnz(II),1);
      z2(II) = z2(II) + 1e-5*rand(nnz(II),1);
    end

    cx(I) = x2;
    cy(I) = y2;
    cz(I) = z2;
    all_dists(I) = dists;
    all_distp(I) = distp;
  end

  cpx(bdy) = cx;
  cpy(bdy) = cy;
  cpz(bdy) = cz;

  dist = sqrt((cpx-xx).^2 + (cpy-yy).^2 + (cpz-zz).^2);
end
