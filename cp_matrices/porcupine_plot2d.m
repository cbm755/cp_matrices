function porcupine_plot(x, y, cpx, cpy, fignum)

  if nargin < 5
    fignum = figure();
  end

  figure(fignum); clf; hold on;
  xlabel('x'); ylabel('y');

  for i=1:length(cpx(:))
    plot([x(i) cpx(i)], [y(i) cpy(i)], 'k-');
    plot(cpx(i),cpy(i), 'r*');
  end

  axis equal