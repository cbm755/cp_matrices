addpath('../cp_matrices');

x = (-2:.1:2)';
y = (-1.5:.1:1.5)';
[xx,yy] = meshgrid(x,y);

radii = [1.4 .2 .8 .4];
[cpx, cpy, dist] = cpRoundedSquare(xx, yy, radii);
figure(1); clf;
porcupine_plot2d(xx, yy, cpx, cpy, 1);
axis equal; axis tight;
plot(xp, yp, 'r-');

figure(2); clf;
pcolor(xx, yy, dist);
hold on;

[xp,yp] = paramRoundedSquare(200, radii);
%figure(3); clf;
plot(xp, yp, 'k-', 'linewidth', 2);
axis equal; axis tight;

