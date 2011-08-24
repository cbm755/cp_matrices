%% Heat equation on a circle
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

% This example solves the heat equation on a 2D circle, with initial
% conditions u = cos(theta), and exact solution u(t) =
% exp(-t)*cos(theta)


%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');
%addpath('/home/wen/closestpoint/cp_matrices');

% add functions for finding the closest points
addpath('../cp_matrices');
%addpath('/home/wen/closestpoint/surfaces');


%%
% 2D example on a circle
% Construct a grid in the embedding space

dx = 0.1;                   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
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
cpxg = cpx(:); cpyg = cpy(:);


%% Banding: do calculation in a narrow band around the sphere
dim = 2;  % dimension
p = 3;    % interpolation order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(dist <= bw*dx);

% store points not in band, for plotting
notband = setdiff(1:nx*ny, band);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band);



%% Function u in the embedding space
% u is a function defined on the grid (eg heat if solving the heat
% equation)

u = zeros(nx,ny);
% assign some initial value (using initial value of cos theta)
[th, r] = cart2pol(xx,yy);
u = cos(th);

% this makes u into a vector, containing only points in the band
u = u(band);

initialu = u;       % store initial value


%% Plot

figure(1);
% make a full matrix for plotting
uplot = zeros(nx,ny);
uplot(band) = u;
uplot(notband) = NaN;

% plot
pcolor(x1d,y1d,uplot);

% make plot look pretty
caxis([-1.1 1.1]);
axis equal; colorbar;
title('initial value on embedding space');
xlabel('x'); ylabel('y');


%% Construct an interpolation matrix for closest point

% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

E = interp2_matrix_band(x1d, y1d, cpxg, cpyg, p, band);

% e.g., closest point extension:
%u = E*u;

%% Create Laplacian matrix for heat equation

order = 2;  % Laplacian order: bw will need to increase if changed
L = laplacian_2d_matrix(x1d,y1d, order, band, band);


%% Construct an interpolation matrix for plotting on circle

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,100);
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix_band(x1d, y1d, xp, yp, p, band);

figure(2); set(gcf,'Position', [410 700 800 800]);


%% Time-stepping for the heat equation

Tf = 1;
dt = 0.2*dx^2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

for kt = 1:numtimesteps
    % explicit Euler timestepping
    unew = u + dt*L*u;

    % closest point extension
    u = E*unew;

    t = kt*dt;
    % plot over computation band
    if ( (kt < 10) | (mod(kt,10) == 0) | (kt == numtimesteps) )
        figure(2);
        subplot(2,1,1); hold off;
        uplot(band) = u;
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
