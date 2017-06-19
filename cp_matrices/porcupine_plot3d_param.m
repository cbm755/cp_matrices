function porcupine_plot3d_param(x, y, z, cpx, cpy, cpz, bdy, paramf, fignum)
%PORCUPINE_PLOT3D_PARAM  Closest Point visualization
%
%   figure(fignum);
%   porcupine_plot3d_param(x,y,z, cpx,cpy,cpz, bdy, paramf, fignum);
%      paramf is the name of function that return a parameterization
%      fignum must already exist.  Works for codim-1 (surfaces) and
%      codim-2 (curves/filaments).
%
%   pass empty to paramf to skip drawing the surface.

  onlyPlotSome = false;

  if (isempty(bdy))
    bdy = zeros(size(cpx));
  end
  figure(fignum); clf;
  if (~ isempty(paramf))
    [xp,yp,zp] = paramf(64);
    if isvector(xp)
      plot3(xp,yp,zp, 'g-', 'linewidth',4);
    else
      surf(xp,yp,zp,zeros(size(xp)));
      shading flat
    end
  end
  hold on;
  xlabel('x'); ylabel('y'); zlabel('z');
  axis equal
  bdycolors = {'k' 'r' 'b' 'm' 'c' 'y'};
  uniqbdy = unique(bdy);
  for c = 1:length(uniqbdy)
    I = (bdy == uniqbdy(c));
    if (onlyPlotSome)
      I = find(I);
      n = numel(I);
      % random
      r = randperm(n); r = r(1:min(500,n));
      % or take the first 500
      %r = 1:min(500, n);
      I = I(r);
    end
    plot3(cpx(I), cpy(I), cpz(I), [bdycolors{c} '.']);
    plot3([x(I) cpx(I)]', [y(I) cpy(I)]', [z(I) cpz(I)]', [bdycolors{c} '-']);
  end

  axis equal
  grid on
