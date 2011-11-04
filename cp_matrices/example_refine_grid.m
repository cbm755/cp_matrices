%% ** Grid Refinement Demo **
% Given a banded coarse grid, we can refine it to finer and finer
% grids MUCH faster than generating from scratch and without building
% big meshgrids (the only fully nD embedding space grid will be at the
% coarest level---which can be negliable).

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');

% add functions for finding the closest points
addpath('../surfaces');

tic

dx = 0.2;  % grid size
% make vectors of x, y, positions of the grid
y1d = (-2:dx:2)';
x1d = ((-2-1*dx):dx:(2+1*dx))';
%x1d = -1.8:dx:1.8;
%y1d = x1d;


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
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band (can discard others)
cpxg = cpx(band); cpyg = cpy(band);
xg = xx(band); yg = yy(band);
distg = dist(band);

time_ini = toc

%% Refine the grid
M = 4;  % how many times to refine
[band2,xg2,yg2,cpx2,cpy2,dist2,dx2,x1d2,y1d2] = ...
    refine_grid(M, cpf, dx, x1d, y1d, bw, band, distg);


%% An alternative approach
% here we keep all the levels of refinement by doing the looping
% manually.
NLevels = 5;

% first level is what we already computed
k = 1;
a_band{k} = band;
a_xg{k} = xg;       a_yg{k} = yg;
a_cpx{k} = cpxg;    a_cpy{k} = cpyg;
a_x1d{k} = x1d;     a_y1d{k} = y1d;
a_dist{k} = distg;
a_dx{k} = dx;

% loop to refine
for k=2:NLevels

  %a_dx{k} = a_dx{k-1}/2;
  [a_band{k}, a_xg{k}, a_yg{k}, a_cpx{k}, a_cpy{k}, a_dist{k}, a_dx{k}, a_x1d{k}, a_y1d{k}] = ...
      refine_grid(1, cpf, a_dx{k-1}, a_x1d{k-1}, a_y1d{k-1}, bw, a_band{k-1}, a_dist{k-1});
end



%% make some plots
bnds = [-2 2 -1.5 1.5];

[xp, yp] = paramf(512);
for k=1:NLevels
  figure(k); clf;
  u = sin(a_xg{k});
  plot2d_compdomain(u,a_xg{k},a_yg{k}, a_dx{k}, a_dx{k}, k)
  axis(bnds);
  plot(xp, yp, 'k-', 'linewidth', 2);
  title(['grid level ' num2str(k) ', dx = ' num2str(a_dx{k})]);
  xlabel('x'); ylabel('y');
end
