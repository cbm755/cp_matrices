function porcupine_plot3d_param(x, y, z, cpx, cpy, cpz, bdy, paramf, fignum)
%PORCUPINE_PLOT3D_PARAM  Closest Point visualization
%
%   figure(fignum);
%   porcupine_plot3d_param(x,y,z, cpx,cpy,cpz, bdy, paramf, fignum);
%      paramf is the name of function that return a parameterization
%      fignum must already exist.  Works for codim-1 (surfaces) and
%      codim-2 (curves/filaments).

  [xp,yp,zp] = paramf(64);
  figure(fignum); clf;
  if isvector(xp)
    plot3(xp,yp,zp, 'g-', 'linewidth',4);
  else
    surf(xp,yp,zp,zeros(size(xp)));
    shading flat
  end
  hold on;
  xlabel('x'); ylabel('y'); zlabel('z');
  N = length(x(:));
  %N = 1000;
  bdycolorset = {'r' 'b' 'm' 'c' 'y'};
  for j=1:N
    %i = ceil(rand*NN);
    i = j;
    %if (abs(z(i)) < 0.2)
      if (bdy(i) == 0)
        plot3([x(i) cpx(i)], [y(i) cpy(i)], [z(i) cpz(i)], 'k-');
        plot3(cpx(i),cpy(i),cpz(i), 'k.');
      else
        h1 = plot3([x(i) cpx(i)], [y(i) cpy(i)], [z(i) cpz(i)], '-');
        h2 = plot3(cpx(i),cpy(i),cpz(i), '-');
        set(h1, 'color', bdycolorset{bdy(i)});
        set(h2, 'color', bdycolorset{bdy(i)});
      end
      %end
  end

  axis equal
  grid on
