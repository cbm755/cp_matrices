function [err_inf, err_L2] = poisson_BeanCurveInterior_Robin(dx)

%% Poisson Equation in a unit disk: $\Delta u - \epsilon*u= f$, 
% Dirichlet B.C.s: $ u = g_Dirichlet $

if nargin < 1
    dx = 0.1;
end
x1d = (-1-3*dx:dx:1+3*dx);
y1d = x1d;

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy] = meshgrid(x1d,y1d);
xx = xx(:); yy = yy(:);

%% Find closest points on the surface
cpf = @(x,y) cpBeanInterior(x,y);

[cpx, cpy, dist, bdy, ~, nx, ny] = cpf(xx, yy);

cpx_bar = cpx; cpy_bar = cpy; dist_bar = dist;
cpx_bar(bdy) = 2*cpx(bdy)-xx(bdy); cpy_bar(bdy) = 2*cpy(bdy)-yy(bdy);
dist_bar(bdy) = 2*dist_bar(bdy);

%% Banding: do calculation in a narrow band around the circle
dim = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) < bw*dx);

% store closest pionts in the band;
cpx = cpx(band); cpy = cpy(band); dist = dist(band);
cpx_bar = cpx_bar(band); cpy_bar = cpy_bar(band); dist_bar = dist_bar(band);
xg = xx(band); yg = yy(band);
nx = nx(band); ny = ny(band);
bdy = bdy(band);

%% Construct interpolation matrices for closest point
% This creates a matrices which interpolate data from the grid x1d y1d, onto
% the points cpx cpy, or cpx_bar cpy_bar.
disp('Constructing interpolation and laplacian matrices');


%% Create Laplacian matrix for Poisson equation
order = 2;
L = laplacian_2d_matrix(x1d,y1d, order, band, band);

%% Exact solution and right hand side function

epsilon = 0;
q = 6; k = 2;
uexactfn = @(x,y) x.^q.*y.^q + x.^q + y.^q + sin(k*pi*x) + sin(k*pi*y) + cos(k*pi*x) + cos(k*pi*y);
rhsfn = @(x,y) q*(q-1)*x.^(q-2).*y.^(q-2).*(x.^2+y.^2) + q*(q-1)*(x.^(q-2) + y.^(q-2)) - k^2*pi^2*(uexactfn(x,y)-x.^q-y.^q-x.^q.*y.^q) - epsilon*uexactfn(x,y);
uxfn = @(x,y) q*x.^(q-1).*y.^q + q*x.^(q-1) + k*pi*cos(k*pi*x) - k*pi*sin(k*pi*x);
uyfn = @(x,y) q*x.^q.*y.^(q-1) + q*y.^(q-1) + k*pi*cos(k*pi*y) - k*pi*sin(k*pi*y);
gammafn = @(x,y) 2 + x + y;

%% Build the coeffecient matrix for Poisson equation consisting of two parts:

lbdy = logical(bdy);
I = speye(size(L));
M = ( L - epsilon*I );

rhs = rhsfn(xg,yg);

E = interp2_matrix(x1d,y1d,cpx(bdy),cpy(bdy),p,band);
Ebar = interp2_matrix(x1d,y1d,cpx_bar(bdy),cpy_bar(bdy),p,band);
cpx_double = 2*cpx_bar(bdy) - cpx(bdy);
cpy_double = 2*cpy_bar(bdy) - cpy(bdy);
Edouble = interp2_matrix(x1d,y1d,cpx_double,cpy_double,p,band);
cpx_triple = 2*cpx_double - cpx_bar(bdy);
cpy_triple = 2*cpy_double - cpy_bar(bdy);
Etriple = interp2_matrix(x1d,y1d,cpx_triple,cpy_triple,p,band);
cpx_quadruple = 2*cpx_triple - cpx_double;
cpy_quadruple = 2*cpy_triple - cpy_double;
Equadruple = interp2_matrix(x1d,y1d,cpx_quadruple,cpy_quadruple,p,band);
cpx_quintuple = 2*cpx_quadruple - cpx_triple;
cpy_quintuple = 2*cpy_quadruple - cpy_triple;
Equintuple = interp2_matrix(x1d,y1d,cpx_quintuple,cpy_quintuple,p,band);


% Robin:
ng = nnz(lbdy);
D = spdiags(dist(lbdy),0,ng,ng);
GAMMA = spdiags(gammafn(cpx(lbdy),cpy(lbdy)),0,ng,ng);

% 1st order approximation to the b.c.
% M(lbdy,:) = D*I(lbdy,:)/dx^2 + GAMMA*(-E + I(lbdy,:))/dx^2;

% 2nd order
 M(lbdy,:) = D*(I(lbdy,:) + Ebar)/2/dx^2 + GAMMA*(-Ebar/2 + I(lbdy,:)/2)/dx^2;

% test
% M(lbdy,:) = D*(I(lbdy,:) + Ebar)/2/dx^2 + + GAMMA*(Edouble/6 - Ebar + E/2 + I(bdy,:)/3) / dx^2;

% 3rd order
% M(lbdy,:) = D*(I(lbdy,:) + 3*Ebar - Edouble) / (3*dx^2) + GAMMA*(Edouble/6 - Ebar + E/2 + I(bdy,:)/3) / dx^2;

% quartic interp, 4th order
% M(lbdy,:) = (-Etriple/12 + Edouble/2 - 1.5*Ebar + 5*E/6 + I(bdy,:)/4) /dx^2;

 
 rhs(lbdy) = ( uexactfn(cpx(lbdy),cpy(lbdy)) +  ...
     gammafn(cpx(lbdy),cpy(lbdy)).*( nx(lbdy).*uxfn(cpx(lbdy),cpy(lbdy)) + ny(lbdy).*uyfn(cpx(lbdy),cpy(lbdy)) ) ).*dist(lbdy)/dx^2;

%  lambda = 2*dim/dx^2;
%  M(lbdy,:) = M(lbdy,:) + lambda*(I(bdy,:) + Ebar) / 2;
%  rhs(lbdy) = rhs(lbdy) + lambda*g_Dirichlet(theta(lbdy));
 
 %% Set up rhs with terms arising from BCs.
tic;
u = M \ rhs;
poisson_time = toc;

%% Compare with exact solution:
error = u - uexactfn(xg,yg);
err_inf = norm(error(~bdy),inf);
err_L2 = norm(error(~bdy),2)*dx;

uexact = uexactfn(xg,yg);
Trun = M*uexact - rhs;
T1 = norm(Trun(~bdy),inf);
T2 = norm(Trun,inf);
[T1,T2]
end