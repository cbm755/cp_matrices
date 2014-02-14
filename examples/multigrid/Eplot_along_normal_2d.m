function [E] = Eplot_along_normal_2d(x1d, y1d, band, x,y,nx,ny,dx,N)
%% build extension matrix for plotting along a certain normal
% x1d, y1d, band: standard input needed for building E matrix
% (x,y): point on the surface
% (nx,ny): unit normal at (x,y)
% dx: meshsize of the embedding grids
% N: At how many points we want to plot.

bw = 3;
x1 = x - bw*dx*nx;
x2 = x + bw*dx*nx;
y1 = y - bw*dx*ny;
y2 = y + bw*dx*ny;

xx = linspace(x1,x2,N)';
yy = linspace(y1,y2,N)';

E = interp2_matrix(x1d,y1d,xx,yy,3);
E = E(:,band);

end
