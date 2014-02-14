%% Heat equation on a circle
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

% This example demonstrates two bands as in the implicit CP paper
% [Macdonald, Ruuth 2009]


%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');


%%
% 2D example on a circle
% Construct a grid in the embedding space

dx = 0.05;   % grid size

% make vectors of x, y, positions of the grid
% x1d = (-2:dx:2)';
% y1d = x1d;
% z1d = x1d;



%% Find closest points on the surface

% meshgrid is only needed for finding the closest points, not afterwards
% [xx yy zz] = meshgrid(x1d, y1d, z1d);
% function cpCircle for finding the closest points on a circle
% [cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);

% using refine_grid to save memory
dx_c = 0.2;
x1d_c = (-2:dx_c:2)';
y1d_c = x1d_c;
z1d_c = x1d_c;

[xx_c yy_c zz_c] = meshgrid(x1d_c,y1d_c,z1d_c);
[cpx_c, cpy_c, cpz_c, dist_c] = cpSphere(xx_c, yy_c, zz_c);

%% Banding: do calculation in narrow bands
% We start by defining an "initial band" which we'll refine later.
% Typically, in practice, we have some scattered input which we'll
% assume satisfies the bandwidth formula below (e.g., output from
% tri2cp)
dim = 3;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed
fd_stenrad = order/2;  % Finite difference stencil radius
% The formula for bw is found in [Ruuth & Merriman 2008] and the
% 1.0002 is a safety factor.
bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((fd_stenrad+(p+1)/2)^2));
% start with a rough band and refine later, here just find the
% indicies of all points within bandwidth of the surface.
band_init_c = find(abs(dist_c) <= bw*dx_c);
distg_c = dist_c(band_init_c);

M = round(log(dx_c/dx)/log(2));
[band_init, xinit, yinit, zinit, cpxinit, cpyinit, cpzinit, distg, bdyg, dx, x1d, y1d, z1d] = ...
refine_grid(M, @cpSphere, dx_c, x1d_c, y1d_c, z1d_c, bw, band_init_c, distg_c);

% the corresponding closest points
% cpxinit = cpx(band_init); cpyinit = cpy(band_init); cpzinit = cpz(band_init);
% xinit = xx(band_init); yinit = yy(band_init); zinit = zz(band_init);

[L,E,R, innerband,outerband, innerbandfull,outerbandfull, innerInOuter] = ...
ops_and_bands3d_test(x1d,y1d,z1d,xinit,yinit,zinit,cpxinit,cpyinit,cpzinit,band_init,p,order);

disp('letting E only do extension to points in oband\iband')
tic;
E(innerInOuter,:) = speye(length(innerband));  %  R(:,innerInOuter) = I.
toc;
disp('done')


% closest points of the inner band
cpxgin = cpxinit(innerband);  
cpygin = cpyinit(innerband);  
cpzgin = cpzinit(innerband);  

%% Function u in the embedding space, initial conditions
% u is a function defined on the grid (e.g. heat)
[th, phi, r] = cart2sph(cpxgin, cpygin, cpzgin);
u0 = cos(phi+pi/2);

%% Diagonal splitting for iCPM
M = lapsharp_unordered(L, E, R);

% plotting grid on sphere, based on parameterization
[xp,yp,zp] = sphere(64);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1, yp1, zp1);

% Eplot interpolate data onto the plotting grid
Eplot = interp3_matrix_test( x1d, y1d, z1d, xp1, yp1, zp1, p);
Eplot = Eplot(:,innerbandfull);

% after building matrices, don't need this set
% ICPM2009BANDINGCHECKS = 0;



% for making full dense plots
%ufull = zeros(ny,nx);
%[th, r] = cart2pol(xx,yy);
%ufull = cos(th);


%% Time-stepping for the heat equation

Tf = 1;
dt = dx/10;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps
err_plot = zeros(numtimesteps,1);

I = speye(size(M));
A = I - dt*M;

uexactfn = @(t,phi) exp(-2*t)*cos(phi+pi/2);

figure(1);
u = u0;
tic
for kt = 1:numtimesteps
  %% implicit backward euler
  [unew flag] = gmres(A,u,[],1e-6);
  u = unew;

  t = kt*dt;
  sphplot = Eplot*u;
  error_sphere_inf = max(abs( uexactfn(t,phi_plot) - sphplot )) / max(abs(uexactfn(t,phi_plot)));
  err_plot(kt) = error_sphere_inf;
  %plotting
   if ( (kt < 5) || (mod(kt,200) == 0) || (kt == numtimesteps) )
       [t dt dx error_sphere_inf]
%       % plot value on sphere
%       figure(2); clf;
%       sphplot = reshape(sphplot, size(xp));
%       surf(xp, yp, zp, sphplot);
%       title( ['soln at time ' num2str(t) ', kt= ' num2str(kt)] );
%       xlabel('x'); ylabel('y'); zlabel('z');
%       axis equal; shading interp;
%       colorbar;
%       pause(0.001);
   end

end
t_two_bands = toc
