function porcupine_plot3d_param(x, y, z, cpx, cpy, cpz, bdy, paramf, fignum)
%PORCUPINE_PLOT3D_PARAM  Closest Point visualization
%
%   figure(fignum);
%   porcupine_plot3d_param(x,y,z, cpx,cpy,cpz, bdy, paramf, fignum);
%      paramf is the name of function that return a parameterization
%      fignum must already exist


%if nargin < 8
%    fignum = figure();
%  end
%  if nargin < 9
%    N = length(cpx(:));
% end

  figure(fignum); clf; hold on;
  xlabel('x'); ylabel('y'); zlabel('z');
  [xp,yp,zp] = paramf(64);
  surf(xp,yp,zp,zeros(size(xp)));
  shading flat
  %xlabel('x'); ylabel('z');

  N = length(x(:));
  %N = 1000;
  for j=1:N
    %i = ceil(rand*NN);
    i = j;
    %if (abs(y(i)) < 0.1)
      if (bdy(i) == 0)
        plot3([x(i) cpx(i)], [y(i) cpy(i)], [z(i) cpz(i)], 'k-');
        plot3(cpx(i),cpy(i),cpz(i), 'k.');
      else
        plot3([x(i) cpx(i)], [y(i) cpy(i)], [z(i) cpz(i)], 'r-');
        plot3(cpx(i),cpy(i),cpz(i), 'r');
      end
      %plot([x(i) cpx(i)], [z(i) cpz(i)], 'k-');
      %plot(cpx(i),cpz(i), 'r*');
      %bend
  end

  axis equal