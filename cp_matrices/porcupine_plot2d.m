function porcupine_plot2d(x, y, cpx, cpy, bdy, paramf, fignum)
%PORCUPINE_PLOT2D  A plot of closest points
%   porcupine_plot2d(x, y, cpx, cpy)
%
%   porcupine_plot2d(x, y, cpx, cpy, bdy)
%      Colour the closest points according to "bdy", an array
%      the same size as cpx, which is nonzero for points whose
%      cp is on the boundary.
%
%   porcupine_plot2d(x, y, cpx, cpy, bdy, paramf)
%      paramf is a function handle that returns a parameterization.
%      Pass empty skip drawing the surface.
%
%   porcupine_plot2d(x, y, cpx, cpy, [], [], fignum)
%      Use figure(fignum); it will be cleared.
%      Pass empty to use the current figure.
%

  if (nargin < 5)
    bdy = [];
  end
  if (nargin < 6)
    paramf = [];
  end
  if (nargin < 7)
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
    [xp, yp] = paramf(256);
    plot(xp,yp, 'g-', 'linewidth', 4);
  end
  hold on;
  xlabel('x'); ylabel('y');

  for i=1:length(cpx(:))
    if (bdy(i) == 0)
      col = 'k';
    elseif (bdy(i) == 1)
      col = 'r';
    elseif (bdy(i) == 2)
      col = 'b';
    elseif (bdy(i) == 3)
      col = 'g';
    elseif (bdy(i) == 4)
      col = 'c';
    elseif (bdy(i) == 5)
      col = 'm';
    else
      col = 'y';
    end
    plot([x(i) cpx(i)], [y(i) cpy(i)], [col '-']);
    plot(cpx(i),cpy(i), [col '*']);
  end

  axis equal

