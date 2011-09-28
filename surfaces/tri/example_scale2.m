% this is used to scale the shape before running the C script
addpath('../readply')

[F,V] = plyread('pulley-orig.ply','tri');

% shift and scale
X = (V(:,1) + 5.5)/ 32;
Y = (V(:,2) + 3.5) / 32;
Z = (V(:,3) + 965) / 32;

Vnew = [X Y Z];

figure(1); clf;
trisurf(F, X, Y, Z, sin(8*Z))
shading interp
pause(0)
xlabel('x'); ylabel('y'); zlabel('z');
axis equal



% F-1 for indexing bug?
plywritetri(F-1,Vnew,'pulley_scaled.ply','ascii');
