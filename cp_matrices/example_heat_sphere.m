%% Heat equation on a sphere
% using cp_matrices

%%
% cp_matrices is a folder of useful functions to make implementing the closest point
% method easier. These include closest point extension matrices, and
% differentiation matrices. 

% This example solves the heat equation on the surface of a sphere, with initial
% conditions u = cos(theta)


%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');

% add function for finding the closest points
addpath('../surfaces');


%%
% 3D example on a sphere
% Construct a grid in the embedding space

dx = 0.1;                   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

nx = length(x1d);
ny = length(y1d);
nz = length(z1d);



%% Find closest points on the surface
% For each point (x,y,z), we store the closest point on the sphere
% (cpx,cpy,cpz)

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy zz] = meshgrid(x1d, y1d, z1d);
% function cpSphere for finding the closest points on a sphere
[cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
% make into vectors
cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);

% do calculation only in a band around the sphere
% band is a vector of the indices of the points in the computation band
band = find(dist <= 4.3*dx);
              
% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);



%% Function u in the embedding space
% u is a function defined on the grid (eg heat if solving the heat
% equation)

% assign some initial value (using initial value of cos (8*theta))
[th, phi, r] = cart2sph(xx,yy,zz);
u = cos(4*th);

% this makes u into a vector, containing only points in the band
u = u(band);  

initialu = u;       % store initial value


%% Construct an interpolation matrix for closest point

% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

p = 3;              % interpolation order
E = interp3_matrix_band(x1d, y1d, z1d, cpxg, cpyg, cpzg, p, band);

% closest point extension
u = E*u;

%% Create Laplacian matrix for heat equation

order = 2;         % Laplacian order
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);


%% Construct an interpolation matrix for plotting on sphere

% plotting grid on sphere
[xp,yp,zp] = sphere(64);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);

Eplot = interp3_matrix_band(x1d, y1d, z1d, xp1, yp1, zp1, p, band);

figure(2); set(gcf,'Position', [410 700 800 800]);


%% Time-stepping for the heat equation

dt = 0.25*dx^2;
timesteps = 20;

for t = 1:timesteps

    % explicit Euler timestepping
    unew = u + dt*L*u;
    
    % closest point extension
    u = E*unew;

    % plot value on sphere
    figure(2);
    sphplot = Eplot*u;
    sphplot = reshape(sphplot, size(xp));
    surf(xp, yp, zp, sphplot);
    title('updated value on sphere');
    xlabel('x'); ylabel('y'); zlabel('z');
    caxis([-1.1 1.1]);
    axis equal; shading interp;
    camlight left; colorbar;
    
end;
