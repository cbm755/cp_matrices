%% Poisson Equation in a disk: $\Delta u - \epsilon*u= f$, 
% with Neumann B.C.s: $ \frac{\partial u}{\partial n} = g $

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices')

% add functions for finding the closest points
addpath('../../surfaces')
addpath('../../surfaces/cpPolygon')

dx = 0.1;
x1d = (-2:dx:2);
y1d = x1d;

nx = length(x1d);
ny = length(y1d);

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy] = meshgrid(x1d,y1d);

%% Find closest points on the surface
cpf = @(x,y) cpCircleInterior(x,y);

[cpx, cpy, dist] = cpf(xx, yy);

% cpbar for BCs
[cpx_bar, cpy_bar, dist_bar, bdy] = cpbar_2d(xx, yy, cpf);

%% Banding: do calculation in a narrow band around the circle
dim = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) < bw*dx);

% store closest pionts in the band;
cpx = cpx(band); cpy = cpy(band); dist = dist(band);
cpx_bar = cpx_bar(band); cpy_bar = cpy_bar(band); dist_bar = dist_bar(band);
x = xx(band); y = yy(band);
bdy = bdy(band);

%% Construct interpolation matrices for closest point
% This creates a matrices which interpolate data from the grid x1d y1d, onto
% the points cpx cpy, or cpx_bar cpy_bar.
disp('Constructing interpolation and laplacian matrices');
E = interp2_matrix_test(x1d,y1d, cpx, cpy, p);
E = E(:,band);
E_bar = interp2_matrix_test(x1d,y1d, cpx_bar, cpy_bar, p);
E_bar = E_bar(:,band);

%% Create Laplacian matrix for Poisson equation
order = 2;
L = laplacian_2d_matrix_test(x1d,y1d, order, band, band);

%% Exact solution and right hand side function
k = 1;
epsilon = 1;
uexactfn = @(theta,r) r.^2.*sin(k*theta);
rhsfn = @(theta,r) (4-k^2)*sin(k*theta) - epsilon*uexactfn(theta,r);
g = @(theta) 2*sin(k*theta);
g_Dirichlet = @(theta) sin(k*theta);

%% Build the coeffecient matrix for Poisson equation consisting of two parts:

lbdy = logical(bdy);
lambda = 2*dim/dx^2;
I = speye(size(L));
M = ( L - epsilon*I );

[theta,r] = cart2pol(cpx_bar,cpy_bar);
rhs = rhsfn(theta,r);

% Neumann:
% cp_bar not converge
%M(lbdy,:) = I(lbdy,:) - E_bar(lbdy,:);
%rhs(lbdy) = g(theta(lbdy)).*dist_bar(lbdy);
% cp 1-st order
%M(lbdy,:) = I(lbdy,:) - E(lbdy,:);
%rhs(lbdy) = g(theta(lbdy)).*dist(lbdy);

% Dirichlet:
M(lbdy,:) = (I(lbdy,:) + E_bar(lbdy,:))/2;
rhs(lbdy) = g_Dirichlet(theta(lbdy));

%% Set up rhs with terms arising from BCs.

tic;
u = M \ rhs;
poisson_time = toc;

%% Compare with exact solution:
error = E_bar*u - uexactfn(theta,r);
norm(error(~lbdy),inf)
