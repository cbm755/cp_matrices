function [err_inf,err_L2] = poisson_disk_eg2_modified(dx)

%% Poisson Equation in a unit disk: $\Delta u - \epsilon*u= f$, 
% Neumann B.C.s: $ \frac{\partial u}{\partial n} = g $
% or Dirichlet B.C.s: $ u = g_Dirichlet $

if nargin < 1
    dx = 0.1;
end
x1d = (-2:dx:2);
y1d = x1d;

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy] = meshgrid(x1d,y1d);

%% Find closest points on the surface
cpf = @(x,y) cpCircleInterior(x,y);

[cpx, cpy, dist, bdy] = cpf(xx, yy);

%% Banding: do calculation in a narrow band around the circle
dim = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) < bw*dx);

% store closest pionts in the band;
cpx = cpx(band); cpy = cpy(band); dist = dist(band);
xg = xx(band); yg = yy(band);
bdy = bdy(band);

% find points in the interior of the cicle for extrapolation.
h = dx;
cpx_bar_bdy = cpx(bdy) - h*(xg(bdy)-cpx(bdy))./dist(bdy);
cpy_bar_bdy = cpy(bdy) - h*(yg(bdy)-cpy(bdy))./dist(bdy);

%% Construct interpolation matrices for closest point
% This creates a matrices which interpolate data from the grid x1d y1d, onto
% the points cpx cpy, or cpx_bar cpy_bar.
disp('Constructing interpolation and laplacian matrices');
% E = interp2_matrix(x1d,y1d, cpx, cpy, p);
% E = E(:,band);
% E_bar = interp2_matrix(x1d,y1d, cpx_bar, cpy_bar, p);
% E_bar = E_bar(:,band);

%% Create Laplacian matrix for Poisson equation
order = 2;
L = laplacian_2d_matrix(x1d,y1d, order, band, band);

%% Exact solution and right hand side function

k = 4; q = 10;
epsilon = 1;
uexactfn = @(theta,r) r.^q.*sin(k*theta);
% \Delta_S f = 1/r d/dr(r df/dr) + 1/r^2 d^2 f/d\theta^2
rhsfn = @(theta,r) (q^2-k^2)*r.^(q-2).*sin(k*theta) - epsilon*uexactfn(theta,r);
g = @(theta) q*sin(k*theta);
g_Dirichlet = @(theta) uexactfn(theta,1);


%% Build the coeffecient matrix for Poisson equation consisting of two parts:

lbdy = logical(bdy);
I = speye(size(L));
M = ( L - epsilon*I );

[theta,r] = cart2pol(xg,yg);
rhs = rhsfn(theta,r);


E = interp2_matrix(x1d,y1d,cpx(bdy),cpy(bdy),p,band);
Ebar = interp2_matrix(x1d,y1d,cpx_bar_bdy,cpy_bar_bdy,p,band);
cpx_double = 2*cpx_bar_bdy - cpx(bdy);
cpy_double = 2*cpy_bar_bdy - cpy(bdy);
Edouble = interp2_matrix(x1d,y1d,cpx_double,cpy_double,p,band);
cpx_triple = 2*cpx_double - cpx_bar_bdy;
cpy_triple = 2*cpy_double - cpy_bar_bdy;
Etriple = interp2_matrix(x1d,y1d,cpx_triple,cpy_triple,p,band);
cpx_quadruple = 2*cpx_triple - cpx_double;
cpy_quadruple = 2*cpy_triple - cpy_double;
Equadruple = interp2_matrix(x1d,y1d,cpx_quadruple,cpy_quadruple,p,band);
cpx_quintuple = 2*cpx_quadruple - cpx_triple;
cpy_quintuple = 2*cpy_quadruple - cpy_triple;
Equintuple = interp2_matrix(x1d,y1d,cpx_quintuple,cpy_quintuple,p,band);


% Neumann:

% linear extrapolation, 2nd order
M(lbdy,:) = ( I(lbdy,:) - Ebar ) / dx^2;
 
rhs(lbdy) = ( h+dist(lbdy) ).*g(theta(lbdy))  / dx^2;

 
 %% Set up rhs with terms arising from BCs.
tic;
u = M \ rhs;
poisson_time = toc;

%% Compare with exact solution:
error = u - uexactfn(theta,r);
err_inf = norm(error(~bdy),inf);
err_L2 = norm(error(~bdy),2)*dx;

uexact = uexactfn(theta,r);
Trun = M*uexact - rhs;
T1 = norm(Trun(~bdy),inf);
T2 = norm(Trun,inf);
[T1,T2]
end