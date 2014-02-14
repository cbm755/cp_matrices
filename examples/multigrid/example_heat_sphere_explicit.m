%% Heat equation on a sphere
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.
%
% This example solves the heat equation on the surface of a sphere,
% with initial conditions u = cos(4*theta)


%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');


%%
% 3D example on a sphere
% Construct a grid in the embedding space

dx = 0.025;                   % grid size
dx = 0.2;
dx_c = 0.4;

% make vectors of x, y, positions of the grid
x1d_c = -2.0:dx_c:2.0;
y1d_c = x1d_c;
z1d_c = x1d_c;

nx_c = length(x1d_c);
ny_c = length(y1d_c);
nz_c = length(z1d_c);



%% Find closest points on the surface
% For each point (x,y,z), we store the closest point on the sphere
% (cpx,cpy,cpz)

% meshgrid is only needed for finding the closest points, not afterwards
[xx_c yy_c zz_c] = meshgrid(x1d_c, y1d_c, z1d_c);
% function cpSphere for finding the closest points on a sphere
[cpx_c, cpy_c, cpz_c, dist_c] = cpSphere(xx_c,yy_c,zz_c);
% make into vectors
cpxg_c = cpx_c(:); cpyg_c = cpy_c(:); cpzg_c = cpz_c(:);


%% Banding: do calculation in a narrow band around the sphere
dim = 3;  % dimension
p = 3;    % interpolation order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band_c = find(abs(dist_c) <= bw*dx_c);

% store closest points in the band;
cpxg_c = cpxg_c(band_c); cpyg_c = cpyg_c(band_c); cpzg_c = cpzg_c(band_c);
xg_c = xx_c(band_c); yg_c = yy_c(band_c); zg_c = zz_c(band_c);
distg_c = dist_c(band_c);

M = round(log(dx_c/dx)/log(2));
[band, xg, yg, zg, cpxg, cpyg, cpzg, distg, bdyg, dx, x1d, y1d, z1d] = ...
    refine_grid(M, @cpSphere, dx_c, x1d_c, y1d_c, z1d_c, bw, band_c, distg_c);



%% Function u in the embedding space
% u is a function defined on the grid (eg heat if solving the heat
% equation)

% assign some initial value (using initial value of cos (8*theta))
[th, phi, r] = cart2sph(xg,yg,zg);
u = cos(phi + pi/2);

initialu = u;       % store initial value


%% Construct an interpolation matrix for closest point

% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

E = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, p);
E = E(:,band);

% e.g., closest point extension:
%u = E*u;

%% Create Laplacian matrix for heat equation

order = 2;  % Laplacian order: bw will need to increase if changed
tic
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);
Ltime = toc

%% Construct an interpolation matrix for plotting on sphere

% plotting grid on sphere, based on parameterization
[xp,yp,zp] = sphere(64);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p);
Eplot = Eplot(:,band);

figure(2); set(gcf,'Position', [410 700 800 800]);


%% Time-stepping for the heat equation

Tf = 0.2;
dt = 0.1*dx^2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps
tic
for kt = 1:numtimesteps
    % explicit Euler timestepping
    unew = u + dt*L*u;

    % closest point extension
    u = E*unew;

    t = kt*dt;

    % plot value on sphere
    if (mod(kt,1000) == 0) || (kt < 10) || (kt == numtimesteps)
      figure(2);
      sphplot = Eplot*u;
      
	  err = norm(exp(-2*t)*cos(phi_plot + pi/2)-sphplot,inf) / norm(exp(-2*t)*cos(phi_plot + pi/2),inf);
      [t dt dx err]
      
      sphplot = reshape(sphplot, size(xp));
      surf(xp, yp, zp, sphplot);
      title( ['soln at time ' num2str(t) ', kt= ' num2str(kt)] );
      xlabel('x'); ylabel('y'); zlabel('z');
      %caxis([-1.05 1.05]);   % lock color axis
      axis equal; shading interp;
%      if ~exist OCTAVE_VERSION camlight left;
      colorbar;
      pause(0.001);
%      end
    end
end
t_explicit = toc


E1 = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, 1,band);
E3 = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, 3,band);
I = speye(size(L));
M = E1*L - 6/dx^2*(I-E3);
%dt = 0.125*dx^2;
dt = 1/6*dx^2;
A = I + dt*M;

kk = 1:100;
norm_Ak = zeros(size(kk));
Ak = I;
for i = kk
    Ak1 = A*Ak;
    tmp = (sum(abs(Ak1),2) - sum(abs(Ak),2));
    %dAk = norm(Ak1-Ak,inf);
    dAk = max(max(Ak1)) - max(max(Ak));
    Ak = Ak1;
    norm_Ak(i) = norm(Ak,inf);
    if mod(i,5) == 0
        disp(['i=',num2str(i),'; norm(A^k,inf)=',num2str(norm_Ak(i)), ...
              '; abs row sum diff: ', num2str(max(tmp)),'  ', num2str(min(tmp)), ...
              '; diff: ', num2str(dAk)]);
    end
end