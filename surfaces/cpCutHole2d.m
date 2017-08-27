function [cpx,cpy,dist,bdy] = cpCutHole2d(x,y,cpf,hole_cen,arclen_target,dx,bw,given_cp_data, tol)
%cpCutHole2d  Cut holes in a curve of prescribed arclength
%
%   Suppose "cpf" is closest point function or a handle to one.
%
%     [cpx, cpy, dist, bdy] = cpCutHole2d(x, y, cpf, ...
%         hole_cen, arclen_target, dx, bw, cp_data)
%
%   x, y: coordinates for which we compute the closest points.
%
%   hole_cen = [hx hy] gives the center of the hole.  If this is not on
%   the surface given by cpf, it will be projected onto the surface (by
%   calling cpf).
%
%   Multiple holes are supported: pass each center as a row of "hole_cen"
%   and pass a vector to "arclen_target".
%
%   See "cpCutHole" for important caveats and additional documentation.

  if (nargin < 9)
    tol = 1e-10;
  end

  if (nargin < 8 || isempty(cp_data))
    [cp_data{1:2}] = cpf(x, y);
  end
  cpx = cp_data{1};
  cpy = cp_data{2};

  bdy = zeros(size(cpx));
  [num_holes, temp] = size(hole_cen);
  assert (temp == 2)

  for k=1:num_holes
    xh = hole_cen(k, 1)
    yh = hole_cen(k, 2)
    % project onto surface
    [xh yh] = cpf(xh, yh)

    % look for a sphere radius which gives approx the target arclength
    R_bar = arclen_target(k)/2
    RR = linspace(0.7*R_bar,1.02*R_bar, 100);
    for i = 1:length(RR)
      I = ((cpx-xh).^2 + (cpy-yh).^2 < RR(i)^2);
      assert(nnz(I) > 0, 'empty circle/surface intersection!')
      SASA(i) = sum(I(:))*(dx)/(2*bw);
    end
    %SASA./arclen_target(k)
    [tilde,i] = min(abs(SASA./arclen_target(k)-1))
    hole_rad(k) = RR(i)
    assert (i > 1 && i < length(RR), 'should be in the middle')
    assert (tilde < 0.1)

    if (true)
      %% TODO: temp plotting
      % index of all points whose cp are within the sphere
      R_H = hole_rad(k);
      I = ((cpx-xh).^2 + (cpy-yh).^2 < R_H^2);
      figure(10);
      porcupine_plot2d(x, y, cpx, cpy);
      plot(xh, yh, 'r*', 'linewidth', 5);
      plot(cpx(I), cpy(I), 'r.');
      plot(x(I), y(I), 'bo');
      t = exp(1i*linspace(0,2*pi,64));
      plot(xh + R_bar*real(t), yh + R_bar*imag(t), 'c:');
      plot(xh + R_H*real(t), yh + R_H*imag(t), 'r-');
    end
  end

  [cpx, cpy, dist, bdy] = cpCutHoleRadius2d(x, y, cpf, hole_cen, hole_rad, cp_data, tol)
end
