addpath('../readply')

[Faces, Vertices] = plyread('beacon.ply', 'tri');

xp = Vertices(:,1);
yp = Vertices(:,2);
zp = Vertices(:,3);

up = sin(10*zp).*cos(7*xp).*sin(6*yp);

figure(2); clf;
trisurf(Faces,xp,yp,zp, up);
camlight left

