function porcupine_plot3d(x, y, z, cpx, cpy, cpz, bdy, paramf, fignum)
%PORCUPINE_PLOT3D  Closest Point visualization
%
%   porcupine_plot3d(x,y,z, cpx,cpy,cpz);
%
%   porcupine_plot3d(x,y,z, cpx,cpy,cpz, bdy);
%      Colour the closest points according to "bdy" which is nonzero
%      for points whose cp is on the boundary.  Pass empty to ignore.
%
%   porcupine_plot3d(x,y,z, cpx,cpy,cpz, bdy, paramf);
%      paramf is a function handle that returns a parameterization
%      Works for codim-1 (surfaces) and codim-2 (curves/filaments).
%      Pass empty to skip drawing the surface.
%
%   porcupine_plot3d(..., fignum);
%      Use figure(fignum); it will be cleared.
%      Pass empty to use the current figure.
%

  onlyPlotSome = false;
  onlyPlotSomeRandom = true;

  if (nargin < 7);
    bdy = [];
  end
  if (nargin < 8)
    paramf = [];
  end
  if (nargin < 9);
    fignum = [];
  end

  if (isempty(bdy))
    bdy = zeros(size(cpx));
  end
  if (isempty(fignum))
    fignum = gcf();
  else
    set(0, 'CurrentFigure', fignum);
  end
  clf;
  if (~ isempty(paramf))
    [xp,yp,zp] = paramf(64);
    if isvector(xp)
      plot3(xp,yp,zp, 'g-', 'linewidth',4);
    else
      surf(xp,yp,zp,zeros(size(xp)));
      shading flat
    end
  else
    view(3)   % otherwise "hold on" shows 2D view
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
      if (onlyPlotSomeRandom)
        r = randperm(n); r = r(1:min(500,n));
      else
        or take the first 500
        r = 1:min(500, n);
      end
      I = I(r);
    end
    plot3(cpx(I), cpy(I), cpz(I), [bdycolors{c} '.']);
    plot3([x(I) cpx(I)]', [y(I) cpy(I)]', [z(I) cpz(I)]', [bdycolors{c} '-']);
  end

  axis equal
  grid on
