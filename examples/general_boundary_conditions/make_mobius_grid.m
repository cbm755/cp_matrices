%% GS on mobius strip

clear ; close all ; % clc ; 

rng(1989)


%% Construct a grid in the embedding space

% grid size
dx = 0.1 ;

% make vectors of x, y, positions of the grid
x1d = (-1.8:dx:1.8)' ;
y1d = x1d ;
z1d = x1d ;

[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;

% Find closest points on the surface
R = 1 ; % radius of center circle
T = 0.35 ; % thickness
cpf = @cpMobiusStrip;
[cpbarx, cpbary, cpbarz, dist, bdy] = cpbar_3d(xx, yy, zz, cpf, R, T) ;
[cpx, cpy, cpz, ~, ~] = cpf(xx, yy, zz) ;


save(['mobius_dx=', num2str(dx), '.mat'], ...
     x1d, y1d, z1d, xx, yy, zz, cpx, cpy, cpz, cpbarx, cpbary, cpbarz, ...
     dx, R, T);

