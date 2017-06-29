
x = (-2:.1:2)';
y = (-1.5:.1:1.5)';
[xx,yy] = meshgrid(x,y);

radii = [1.0 .2 .8 .4];
[cpx, cpy, dist] = cpRoundedSquare(xx, yy, radii);
figure(1);
porcupine_plot2d(xx, yy, cpx, cpy);
axis equal; axis tight;
[xp,yp] = paramRoundedSquare(200, radii);
plot(xp, yp, 'r-', 'linewidth', 2);

figure(2); clf;
pcolor(xx, yy, dist);
hold on;

%figure(3); clf;
plot(xp, yp, 'k-', 'linewidth', 2);
axis equal; axis tight;

