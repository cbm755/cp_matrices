%% ** Grid Refinement Demo **
% Given a banded coarse grid, we can refine it to finer and finer
% grids MUCH faster than generating from scratch and without building
% big meshgrids (the only fully nD embedding space grid will be at the
% coarest level---which can be made almost negliable).
%
% TODO: doesn't deal with open surfaces as well as it should: bdy
% not computed.

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');

% add functions for finding the closest points
addpath('../surfaces');


dx = 0.25;  % coarsest grid size
% make vectors of x, y, positions of the grid
y1d = (-2:dx:2)';
x1d = ((-2-1*dx):dx:(2+1*dx))';


%% Find a coarse band of closest points
% meshgrid is only needed at this coarse grid
[xx yy] = meshgrid(x1d, y1d);
%cpf = @cpEllipse;  paramf = @paramEllipse;
%cpf = @cpCircle;  paramf = @paramCircle;
cpf1 = @cpSemicircle;  paramf = @paramSemicircle;  cpf = @(x,y) cpbar_2d(x,y,cpf1);
[cpx, cpy, dist] = cpf(xx,yy);

%% Banding: make a narrow band around the circle
dim = 2;  % dimension
p = 3;    % interpolation order
% "band" is a vector of the indices of the points in the computation
% band.
bw = rm_bandwidth(dim, p);
band = find(abs(dist) <= bw*dx);

% store closest points in the band (can discard others)
cpxg = cpx(band); cpyg = cpy(band);
xg = xx(band); yg = yy(band);
distg = dist(band);

gc = [];
gc.dim = 2;
gc.dx = dx;
gc.x1d = x1d;
gc.y1d = y1d;
gc.cpfun = cpf;
gc.band = band;
gc.x = xg;
gc.y = yg;
gc.cpx = cpxg;
gc.cpy = cpyg;
gc.dist = distg;

%% Refine the grid

NLevels = 5;
g = {};
g{1} = gc;
for i=2:NLevels
  g{i} = refine_cpgrid_bw(g{i-1}, bw);
end


%% make some plots
bnds = [-2 2 -1 2];

[xp, yp] = paramf(512);
for k=1:NLevels
  figure(k); clf;
  u = sin(g{k}.x);
  plot2d_compdomain(u, g{k}.x, g{k}.y, g{k}.dx, g{k}.dx, k);
  axis(bnds);
  plot(xp, yp, 'k-', 'linewidth', 2);
  title(['grid level ' num2str(k) ', dx = ' num2str(g{k}.dx)]);
  xlabel('x'); ylabel('y');
end
