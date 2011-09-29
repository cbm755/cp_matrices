function porcupine_plot2d(x, y, cpx, cpy, fignum, bdy)
%PORCUPINE_PLOT2D  A plot of closest points
%   porcupine_plot2d(x, y, cpx, cpy, fignum, bdy)

  if nargin < 5
    fignum = figure();
  end
  if nargin < 6
    bdy = zeros(size(cpx));
  end

  %figure(fignum);
  set(0, 'CurrentFigure', fignum);
  clf; hold on;
  xlabel('x'); ylabel('y');

  for i=1:length(cpx(:))
    if (bdy(i) == 0)
      col = 'k';
    elseif (bdy(i) == 1)
      col = 'r';
    elseif (bdy(i) == 2)
      col = 'b';
    else
      col = 'y';
    end
    plot([x(i) cpx(i)], [y(i) cpy(i)], [col '-']);
    plot(cpx(i),cpy(i), [col '*']);
  end

  axis equal

