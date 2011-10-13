function poly = initPoly
poly = [];

figure; 
plot([-1 1 1 -1 -1],[-1 -1 1 1 -1],'k');
axis([-1.5 1.5 -1.5 1.5]);
axis equal;

h=impoly;
poly=wait(h);