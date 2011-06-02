% this is used to scale the shape before running the C script
addpath('../readply')

[F,V] = plyread('bumpy_torus_orig.ply','tri');

% scale
V = V / 10;

X = V(:,1);
Y = V(:,2);
Z = V(:,3);

figure(1); clf;
trisurf(F, X, Y, Z, sin(8*Z))
shading interp
pause(0)

% F-1 for indexing bug?
plywritetri(F-1,V,'bumpy_torus_scaled.ply','ascii');
