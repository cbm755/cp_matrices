%% Gray--Scott reaction-diffusion on a triangulated pig.
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

warning('use example_RD_tri.m instead');

% adjust as appropriate
addpath('../../cp_matrices');
addpath('../../surfaces');
addpath('../../surfaces/tri_pig');
addpath('../../surfaces/readply');

loaddata = 1;

% ply file contains the triangles (for plotting, see below)
PlyFile = 'pig_loop2.ply';

dim = 3;
p = 3;  % degree interp
order = 2;  % laplacian order, griddata hardcoded for 2.
% NOTE: griddata needs to be generated for each choice of dx (by the
% tri2cp C code).  its loaded below.  dim, p, order, x1d, y1d, z1d
% must match this code too.

w = 1;
p_f2c = 1;
p_c2f = 1;
MAX = 1;
n1 = 2;
n2 = 2;

dx = 0.025;
dx_coarsest = 0.1;

n_level = round(log(dx_coarsest/dx)/log(2)) + 1;

a_band = cell(n_level,1);
a_xcp = cell(n_level,1);
a_ycp = cell(n_level,1);
a_zcp = cell(n_level,1);
a_x1d = cell(n_level,1);
a_y1d = cell(n_level,1);
a_z1d = cell(n_level,1);
a_xg = cell(n_level,1);
a_yg = cell(n_level,1);
a_zg = cell(n_level,1);

Lc = cell(n_level,1);
Ec = cell(n_level,1);
Mc = cell(n_level,1);

a_dist = cell(n_level,1);
a_bdyg = cell(n_level,1);

dx_tmp = dx;

for m = 1:1:n_level

a_x1d{m} = (-2.0:dx_tmp:2.0)';
a_y1d{m} = a_x1d{m};
a_z1d{m} = a_x1d{m};

nx = length(a_x1d{m});
ny = length(a_y1d{m});
nz = length(a_z1d{m});
% the griddata file assumes something like this for ordering, but
% its not necessary to actually do the meshgrid.
%[x3d,y3d,z3d]=meshgrid(x1d,y1d,z1d);
%[cpx,cpy,cpz,dist] = cpSphere(x3d,y3d,z3d);

% griddata: processed from ply file by tri2cp (C code).  Contains
% only points in the narrow band.
GD = load( ['~/svn/closestpoint/surfaces/tri_pig/pig_loop2_griddata_p' num2str(p) ...
           '_dx' num2str(dx_tmp) '.txt'] );
% plus one b/c they're C indices
i = GD(:,1) + 1;
j = GD(:,2) + 1;
k = GD(:,3) + 1;
a_dist{m} = GD(:,4);
a_xcp{m} = GD(:,5);
a_ycp{m} = GD(:,6);
a_zcp{m} = GD(:,7);
a_xg{m} = GD(:,8);
a_yg{m} = GD(:,9);
a_zg{m} = GD(:,10);


%%sanity checks
xtest = a_x1d{m}(1) + (i-1)*dx_tmp;
ytest = a_y1d{m}(1) + (j-1)*dx_tmp;
ztest = a_z1d{m}(1) + (k-1)*dx_tmp;
if ( (norm(xtest - a_xg{m}, inf) > 1e-14) ||  ...
     (norm(ytest - a_yg{m}, inf) > 1e-14) || ...
     (norm(ztest - a_zg{m}, inf) > 1e-14) )
  error('sanity fail');
end
%[temp,N] = size(GD);
%for c=1:N

%end

% here is one place where meshgrid comes in: note ordering here.
band = sub2ind([ny,nx,nz], j,i,k);
a_band{m} = band;

build_matrices = true;
if (build_matrices)
  %% discrete operators
  disp('building laplacian and interp matrices');
  tic;
  Lc{m} = laplacian_3d_matrix_test(a_x1d{m}, a_y1d{m}, a_z1d{m}, order, band,band);
  toc;
  tic;
  Ec{m} = interp3_matrix_test(a_x1d{m}, a_y1d{m}, a_z1d{m}, a_xcp{m}, a_ycp{m}, a_zcp{m}, p);
  Ec{m} = Ec{m}(:,band);
  toc;
  % iCPM matrix
  Mc{m} = lapsharp(Lc{m},Ec{m});

  
  %Eplot0 = interp3_matrix_band(x1d,y1d,z1d, xp,yp,zp, 0, band);
end

dx_tmp = dx_tmp*2;

end

L = Lc{1};
E = Ec{1};

disp('building transform matrices to do restriction and prolongation ...')
[TMf2c, TMc2f] = helper_set_TM3d(a_x1d, a_y1d, a_z1d, a_xcp, a_ycp, a_zcp, a_band, a_bdyg, p_f2c, p_c2f);
disp('done')

%% plotting grid
[Faces, Vertices] = plyread(PlyFile, 'tri');
xp = Vertices(:,1);
yp = Vertices(:,2);
zp = Vertices(:,3);
Eplot = interp3_matrix_test(a_x1d{1}, a_y1d{1}, a_z1d{1}, xp,yp,zp, 3);
Eplot = Eplot(:,a_band{1});
Eplot1 = interp3_matrix_test(a_x1d{1}, a_y1d{1}, a_z1d{1}, xp,yp,zp, 1);
Eplot1 = Eplot1(:,a_band{1});
  
uexact = zeros(size(xp));

% parameters and functions for Gray--Scott
% 120 works with 0.025
F = 0.054;  k = 0.063;  nuu = 1/120^2;  nuv = nuu/2;
f = @(u,v) (-u.*v.*v  +  F*(1-u));
g = @(u,v) ( u.*v.*v  -  (F+k)*v);
  
%% Do some calculation
% Here its eigenvalues/eigenfunctions but could be timestepping
% or whatever using the cp_matrices.

xg = a_xg{1}; yg = a_yg{1}; zg = a_zg{1};
[th,r,temp] = cart2pol(xg,yg,zg);
%pert = 0.5*exp(-(6*(zg-.1)).^2) + 0.5*rand(size(xg));
pert = 1*exp(-(6*(zg-0.05*cos(6*th))).^2); % + 0.5*rand(size(xg));
u0 = 1 - pert;
v0 = 0 + 0.5*pert;

Tf = 15000;

% explicit dt for dx = 0.025
% dt = 1;

% implicit dt for dx = 0.025
dt = 10;

numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

lambda_u = 2*nuu/3;
lambda_v = 2*nuv/3;

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

implicit = 1;

u = u0;
v = v0;

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


for kt = 1:numtimesteps
    % explicit Euler timestepping
    if (~implicit)
    tic
    unew = u + dt*f(u,v) + dt*nuu*(L*u);
    vnew = v + dt*g(u,v) + dt*nuv*(L*v);
    toc
    
    % closest point extension
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

     [u2] = gmres(Mu{1}, rhs_u, [], 1e-3);
     [v2] = gmres(Mv{1}, rhs_v, [], 1e-3);

    % u2 = gmg_t(Mu, DMu, Lu, Ec, U, Fu, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w, uexact, Eplot, MAX);
    % v2 = gmg_t(Mv, DMv, Lv, Ec, V, Fv, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w, uexact, Eplot, MAX);

    % u2 = helper_vcycle_t(Mu, DMu, Lu, Ec, U, Fu, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w);
    % v2 = helper_vcycle_t(Mv, DMv, Lv, Ec, V, Fv, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w);

    % u2 = gmg(Mu, Lu, Ec, U, Fu, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w, uexact, Eplot, MAX);
    % v2 = gmg(Mv, Lv, Ec, V, Fv, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w, uexact, Eplot, MAX);

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
    
    if ( (mod(kt,4)==0) || (kt<=10) || (kt==numtimesteps) )
    % plot value on sphere
    figure(1); clf;
    uplot = Eplot*u;
    trisurf(Faces,xp,yp,zp, uplot);
    xlabel('x'); ylabel('y'); zlabel('z');
    title( ['u at time ' num2str(t)] );
    axis equal
    shading interp
    camlight left
    colorbar
    pause(0.001);
    end

end
