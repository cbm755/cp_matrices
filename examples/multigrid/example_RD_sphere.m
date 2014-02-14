%% Gray--Scott reaction-diffusion on sphere, explicit time
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.


% adjust as appropriate
addpath('../../cp_matrices');
addpath('../../surfaces');

loaddata = 1;

if (loaddata == 1)
  %%
% 3D example on a sphere
% Construct a grid in the embedding space

dx = 0.05;                   % grid size
dx_coarsest = 0.4;

dt = dx/8;

% parameters and functions for Gray--Scott
% 120 works with 0.025, nu sets the scale of the patterns: here we
% set it related to dx so that finer grids have finer patterns
% F = 0.054;  k = 0.063;  nuu = 1/(3/dx)^2;  nuv = nuu/2;
% f = @(u,v) (-u.*v.*v  +  F*(1-u));
% g = @(u,v) ( u.*v.*v  -  (F+k)*v);
nuu = 1.25/900;
nuv = 10/900;
a = 3; b = 10.2;
f = @(u,v)  (a - (b+1)*u + u.*u.*v);
g = @(u,v)  (b*u - u.*u.* v);

% make vectors of x, y, positions of the grid
x1d_coarsest = (-3.0:dx_coarsest:3.0)';
y1d_coarsest = x1d_coarsest;
z1d_coarsest = x1d_coarsest;

has_boundary = false;

%% Banding: do calculation in a narrow band around the sphere
dim = 3;  % dimension
p = 2;    % interpolation order
order = 2;
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));


MAX = 1;
n1 = 3;
n2 = n1 + round(log(nuv/nuu)/log(2)/2);
p_f2c = 1;
p_c2f = 1;
w = 1;

cpf = @cpSphere;

disp('building closest point grids for each level of V-Cycle ... ')
[a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d] = ...
    build_mg_cpgrid3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);
disp('done');

n_level = length(a_band);

%lambda_u = 0.5*nuu;
%lambda_v = 0.5*nuv;

lambda_u = 2*nuu/3;
lambda_v = 2*nuv/3;

disp('building cp matrices for each level of V-Cycle ...')
%[Lu, Ec, L1_u, M1_u, Mn_u] = build_mg_cpmatrix3d_test(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order, 1, dt, lambda_u);
%[Lv, Ec, L1_v, M1_v, Mn_v] = build_mg_cpmatrix3d_test(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order, 1, dt, lambda_v);
%[Mu, Lu, Ec, L1_u] = build_mg_cpmatrix3d(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order, 1, dt, lambda_u);
%[Mv, Lv, Ec, L1_v] = build_mg_cpmatrix3d(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order, 1, dt, lambda_v);
[Mc, Lc, Ec] = build_mg_cpmatrix3d(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order);
disp('done');

%L = L1_u;
L = Lc{1};
E = Ec{1};

disp('fix the coefficient matrix for the time dependent problem')
tic;
n_level = length(a_band);
Mu = cell(n_level,1);
Mv = cell(n_level,1);
DMu = cell(n_level,1);
DMv = cell(n_level,1);
Lu = Lc;
Lv = Lc;
for i = 1:1:n_level
    Mu{i} = speye(size(Mc{i})) - lambda_u*dt*Mc{i};
    Mv{i} = speye(size(Mc{i})) - lambda_v*dt*Mc{i};
    DMu{i} = diag(Mu{i});
    DMv{i} = diag(Mv{i});
    Lu{i} = speye(size(Lc{i})) - lambda_u*dt*Lc{i};
    Lv{i} = speye(size(Lc{i})) - lambda_v*dt*Lc{i};
end

toc;
disp('done')

disp('building transform matrices to do restriction and prolongation ...')
[TMf2c, TMc2f] = helper_set_TM3d(a_x1d, a_y1d, a_z1d, a_xcp, a_ycp, a_zcp, a_band, a_bdyg, p_f2c, p_c2f);
disp('done')

cpxg = a_xcp{1};
cpyg = a_ycp{1};
cpzg = a_zcp{1};

Radius = 1;

% plotting grid on sphere, based on parameterization
[xp,yp,zp] = paramSphere(64,Radius);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = cell(n_level,1);
for i = 1:1:n_level
Eplot{i} = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, xp1, yp1, zp1, p);
Eplot{i} = Eplot{i}(:, a_band{i});
end

uexact = zeros(size(xp1));

end
    
%% Do some calculation
% Here its eigenvalues/eigenfunctions but could be timestepping
% or whatever using the cp_matrices.

[th,r,temp] = cart2pol(cpxg,cpyg,cpzg);
pert = 1.6*rand(size(cpxg));
%pert = 0.8*exp(-(1*(cos(0*th)+2).*(r)).^2) + 0.7*rand(size(cpxg));
%pert = 0.8*exp(-(3*r.^2)) + 0.7*rand(size(cpxg));
%pert = 0.5*exp(-(6*(cpzg-.2)).^2) + 0.5*rand(size(cpxg));
%pert = 0.9*exp(-(6*(zg-0.25*(cos(6*th))-.5)).^2) + 0.1*rand(size(xg));
%u0 = 1 - pert;
%v0 = 0 + 0.5*pert;
%pert = 0.1*rand(size(cpxg));
%pert = zeros(size(cpxg));
u0 = a + pert;
v0 = b/a + pert;
u = u0;
v = v0;

%u = zeros(size(cpxg));
%v = zeros(size(cpxg));


implicit = 1

if ~implicit
    dt = 0.1*dx^2;  % this is basically 1 b/c nuu is like dx^2
    %dt = dx^2;
end

Tf = 2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

M1 = speye(size(Mc{1})) - dt*nuu*Mc{1};
M2 = speye(size(Mc{1})) - dt*nuv*Mc{1};
rhs_u = u + dt*f(u,v);
rhs_v = v + dt*g(u,v);
[u1 flag] = gmres(M1,rhs_u);
[v1 flag] = gmres(M2,rhs_v);
    
%unew = u + dt*f(u,v) + dt*nuu*(L*u);
%vnew = v + dt*g(u,v) + dt*nuv*(L*v);
%u1 = E*unew;
%v1 = E*vnew;

U = cell(n_level,1);
V = cell(n_level,1);
Fu = cell(n_level,1);
Fv = cell(n_level,1);
for kt = 2:numtimesteps
    if (~implicit)
    %explicit Euler timestepping
    tic
    unew = u + dt*f(u,v) + dt*nuu*(L*u);
    vnew = v + dt*g(u,v) + dt*nuv*(L*v);
    toc
    %closest point extension
    tic
    u = E*unew;
    v = E*vnew;
    toc
    else

    % IMEX
    tic
    for i = 2:1:n_level
        U{i} = zeros(size(a_band{i}));
        V{i} = zeros(size(a_band{i}));
    end
    %rhs_u = u1 + 0.5*dt*( 3*f(u1,v1) - f(u,v) + nuu*L*u1 );
    %rhs_v = v1 + 0.5*dt*( 3*g(u1,v1) - g(u,v) + nuv*L*v1 );
    %rhs_u = u1 + 0.5*dt*( 3*f(u1,v1) - f(u,v) );
    %rhs_v = v1 + 0.5*dt*( 3*g(u1,v1) - g(u,v) );
    
    rhs_u = ( 4*u1 - u + dt*( 4*f(u1,v1) - 2*f(u,v) ) )/3;
    rhs_v = ( 4*v1 - v + dt*( 4*g(u1,v1) - 2*g(u,v) ) )/3;
	Fu{1} = rhs_u;
    Fv{1} = rhs_v;
    U{1} = Fu{1};
    V{1} = Fv{1};
    %U{1} = u1;
    %V{1} = v1;
    % u2 = M1_u \ rhs_u;
    % v2 = M1_v \ rhs_v;
    
    % [u2] = gmres(Mu{1}, rhs_u, 5, 1e-12);
    % [v2] = gmres(Mv{1}, rhs_v, 5, 1e-12);
    
     u2 = gmg_t(Mu, DMu, Lu, Ec, U, Fu, TMf2c, TMc2f, a_band, a_bdyg, n1, n1, 1, w, uexact, Eplot, MAX);
     v2 = gmg_t(Mv, DMv, Lv, Ec, V, Fv, TMf2c, TMc2f, a_band, a_bdyg, n2, n2, 1, w, uexact, Eplot, MAX);
    
    % u2 = helper_vcycle_t(Mu, DMu, Lu, Ec, U, Fu, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w);
    % v2 = helper_vcycle_t(Mv, DMv, Lv, Ec, V, Fv, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w);
    
    % u2 = gmg(Mu, Lu, Ec, U, Fu, TMf2c, TMc2f, a_band, a_bdyg, n1, n1, 1, w, uexact, Eplot, MAX);
    % v2 = gmg(Mv, Lv, Ec, V, Fv, TMf2c, TMc2f, a_band, a_bdyg, n2, n2, 1, w, uexact, Eplot, MAX);
    
    % u2 = helper_vcycle_M(Mu, Lu, Ec, U, Fu, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w);
    % v2 = helper_vcycle_M(Mv, Lv, Ec, V, Fv, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w);
    
    % u2 = helper_vcycle(Mu, Lu, Ec, U, Fu, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w);
    % v2 = helper_vcycle(Mv, Lv, Ec, V, Fv, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w);

    %     res_u = Eplot{1}*(rhs_u - Mu{1}*u2);
    %     res_v = Eplot{1}*(rhs_v - Mv{1}*v2);
    res_u = (rhs_u - Mu{1}*u2);
    res_v = (rhs_v - Mv{1}*v2);
    resu = norm(res_u,inf)
    resv = norm(res_v,inf)
    u = u1;  v = v1;
    u1 = u2; v1 = v2;
    toc
    end
    
    t = kt*dt;
    
    if ( (mod(kt,50)==0) | (kt<=25) | (kt==numtimesteps) )
    % plot value on sphere
    figure(3); clf;
    sphplot = Eplot{1}*u;
    sphplot = reshape(sphplot, size(xp));
    surf(xp, yp, zp, sphplot);
    xlabel('x'); ylabel('y'); zlabel('z');
    title( ['u at time ' num2str(t)] );
    axis equal
    shading interp
    camlight left
    colorbar
    pause(0.01);
    end

end
