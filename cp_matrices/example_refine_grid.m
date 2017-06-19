%% ** Grid Refinement Demo **
% Given a banded coarse grid, we can refine it to finer and finer
% grids MUCH faster than generating from scratch and without building
% big meshgrids (the only fully nD embedding space grid will be at the
% coarest level---which can be made almost negliable).
%
% Note: this does require that you can call a closest point function
% (rather than just having a samples of one).
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
%cpf = @cpBeanCurve;  paramf = @paramBeanCurve;
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
%cpxg = cpx(:); cpyg = cpy(:);
%xg = xx(:); yg = yy(:);
%distg = dist(:);

%% "pack" the grid into a struct
g = [];
g.dim = 2;
g.dx = dx;
g.x1d = x1d;
g.y1d = y1d;
g.cpfun = cpf;
g.band = band;
g.x = xg;
g.y = yg;
g.cpx = cpxg;
g.cpy = cpyg;
g.dist = distg;

%% Refine the grid
g2 = refine_cpgrid_bw(g, bw);
g3 = refine_cpgrid_bw(g2, bw);
g4 = refine_cpgrid_bw(g3, bw);
g5 = refine_cpgrid_bw(g4, bw);

gg = {g g2 g3 g4 g5};

%% or we could use a loop
%NLevels = 5;
%gg = {};
%gg{1} = g;
%for i=2:NLevels
%  gg{i} = refine_cpgrid_bw(gg{i-1}, bw);
%end

%% make some plots
bnds = [-2 2 -1 2];

[xp, yp] = paramf(512);
for k=1:length(gg)
  figure(k); clf;
  u = sin(gg{k}.x);
  plot2d_compdomain(u, gg{k}.x, gg{k}.y, gg{k}.dx, gg{k}.dx, k);
  axis(bnds);
  plot(xp, yp, 'k-', 'linewidth', 2);
  title(['grid level ' num2str(k) ', dx = ' num2str(gg{k}.dx)]);
  xlabel('x'); ylabel('y');
end
