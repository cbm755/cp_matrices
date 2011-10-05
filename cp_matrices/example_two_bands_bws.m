%% Heat equation on a circle
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

% This example demonstrates two bands as in implicit CP paper
% [Macdonald, Ruuth 2009]
%%%%%%%%%%%%%
% Work in progress!  use at your own risk
%%%%%%%%%%%%%


%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');

% add functions for finding the closest points
addpath('../surfaces');



%%
% 2D example on a circle
% Construct a grid in the embedding space

dx = 0.2;   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2.4:dx:2.4)';
y1d = x1d;

nx = length(x1d);
ny = length(y1d);



%% Find closest points on the surface
% For each point (x,y), we store the closest point on the circle
% (cpx,cpy)

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy] = meshgrid(x1d, y1d);
% function cpCircle for finding the closest points on a circle
[cpx, cpy, dist] = cpCircle(xx,yy);
% make into vectors
%cpxg = cpx(:); cpyg = cpy(:);


%% Banding: do calculation in a narrow band around the sphere
dim = 2;  % dimension
p = 3;    % interpolation order
order = 4;  % Laplacian order: bw will need to increase if changed
fd_stenrad = order/2;  % Finite difference stencil radius
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
%bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));

bw1 = 1.0001*sqrt((dim)*((p+1)/2)^2);
bw2 = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((fd_stenrad+(p+1)/2)^2));
band1 = find(abs(dist) <= bw1*dx);
band2 = find(abs(dist) <= bw2*dx);
% todo: had some slick way to find bw1 indices in bw2
%band1in2 = 

nband1 = length(band1)
nband2 = length(band2)

% store points not in band, for plotting
notband = setdiff(1:nx*ny, band1);

% store closest points in the band;
cpxg1 = cpx(band1); cpyg1 = cpy(band1);
cpxg2 = cpx(band2); cpyg2 = cpy(band2);



%% Function u in the embedding space
% u is a function defined on the grid (eg heat if solving the heat
% equation)

u = zeros(nx,ny);
% assign some initial value (using initial value of cos theta)
[th, r] = cart2pol(xx,yy);
u = cos(th);

% this makes u into a vector, containing only points in the band
u = u(band1);


initialu = u;       % store initial value


%% Plot

figure(1);
% make a full matrix for plotting
uplot = zeros(nx,ny);
uplot(band1) = u;
uplot(notband) = -1;

% plot
pcolor(x1d,y1d,uplot);

% make plot look pretty
%caxis([-1.1 1.1]);
axis equal; colorbar;
title('initial value on embedding space');
xlabel('x'); ylabel('y');


%% Construct an interpolation matrix for closest point

% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');
% TODO: should just make one: have some sort of restriction operator
E1 = interp2_matrix_band(x1d, y1d, cpxg1, cpyg1, p, band1);
E = interp2_matrix_band(x1d, y1d, cpxg2, cpyg2, p, band1);


%% Create Laplacian matrix for heat equation


L = laplacian_2d_matrix(x1d,y1d, order, band1, band2);
Ls = laplacian_2d_matrix(x1d,y1d, order, band2, band2);

M = lapsharp(L, E);

%% Construct an interpolation matrix for plotting on circle

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,100);
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix_band(x1d, y1d, xp, yp, p, band1);

figure(2); set(gcf,'Position', [410 700 800 800]);




%% Time-stepping for the heat equation

Tf = 1;
dt = 0.2*dx^2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

for kt = 1:numtimesteps
    %u2 = E*u;
    %u2short = 
    % explicit Euler timestepping
    unew = E1*u + dt*L*(E*u);

    % closest point extension
    %u = E*unew;
    u = unew;

    t = kt*dt;
    % plot over computation band
    if ( (kt < 5) || (mod(kt,20) == 0) || (kt == numtimesteps) )
        figure(2);
        subplot(2,1,1); hold off;
        uplot(band1) = u;
        pcolor(x1d,y1d,uplot);
        caxis([-1.05 1.05]);
        hold on
        plot(xp,yp,'k-', 'linewidth', 2);
        axis equal;  axis tight;
        colorbar;
        title( ['embedded domain: soln at time ' num2str(t) ...
                ', timestep #' num2str(kt)] );
        xlabel('x'); ylabel('y');

        % plot value on circle
        figure(2); subplot(2,1,2); hold off;
        circplot = Eplot*u;
        plot(thetas, circplot);
        title( ['soln at time ' num2str(t) ', on circle'] );
        xlabel('theta'); ylabel('u');
        hold on;
        % plot analytic result
        plot(thetas, exp(-t)*cos(thetas), 'r--');
        legend('explicit Euler', 'exact answer', 'Location', 'SouthEast');
        error_circ_inf = max(abs( exp(-t)*cos(thetas) - circplot' ))
        pause(0);
    end
end
