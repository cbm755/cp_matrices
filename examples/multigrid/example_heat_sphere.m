%% Heat equation on a sphere
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


x0 = -2;
x1 = 2;
y0 = -2;
y1 = 2;
z0 = -2;
z1 = 2;

%%
% 2D example on a circle
% Construct a grid in the embedding space

dx = 0.05; % grid size
dx_coarsest = 0.2;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';
z1d_coarsest = (z0:dx_coarsest:z1)';

dt = dx/10;
lambda = 1;

dim = 3;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

w = 1;
n1 = 1;
n2 = 1;
nc = 3;

MAX = 1;

p_f2c = 2;
p_c2f = 2;

cpf = @cpSphere;


has_boundary = false;
disp('building closest point grids for each level of V-Cycle ... ')
[a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d] = ...
    build_mg_cpgrid3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);
disp('done');

disp('building cp matrices for each level of V-Cycle ...')
[Mc, Lc, Ec, L1] = build_mg_cpmatrix3d(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order, 1, dt, lambda);
disp('done');

% disp('fix the coefficient matrix for the time dependent problem')
% tic;
% n_level = length(a_band);
% M1 = speye(size(M1)) - lambda*dt*M1;
% Mn = speye(size(Mn)) - lambda*dt*Mn;
% toc;
% disp('done')

disp('building transform matrices to do restriction and prolongation ...')
[TMf2c, TMc2f] = helper_set_TM3d(a_x1d, a_y1d, a_z1d, a_xcp, a_ycp, a_zcp, a_band, a_bdyg, p_f2c, p_c2f);
disp('done')

%uexactfn = @(t, th, phi) exp(-2*t)*cos(phi+pi/2);

%uexactfn = @(t, th, phi) exp(-2*t)*( cos(phi+pi/2) + sin(phi+pi/2).*cos(th) );

%uexactfn = @(t, th, phi) exp(-2*t)*( cos(phi+pi/2) + sin(phi+pi/2).*cos(th) ) + ...
%                         exp(-6*t)*(3*cos(phi+pi/2).^2-1 + cos(2*th).*sin(phi+pi/2).^2);

uexactfn = @(t, th, phi) exp(-2*t)*(cos(th).*sin(phi+pi/2)) + exp(-6*t)*(cos(2*th).*sin(phi+pi/2).^2) + ...
                         exp(-12*t)*(cos(3*th).*sin(phi+pi/2).^3) + exp(-20*t)*(cos(4*th).*sin(phi+pi/2).^4) + ...
                         exp(-30*t)*(cos(5*th).*sin(phi+pi/2).^5) + exp(-42*t)*(cos(6*th).*sin(phi+pi/2).^6);

cpxg = a_xcp{1};
cpyg = a_ycp{1};
cpzg = a_zcp{1};
[th, phi, r] = cart2sph(cpxg,cpyg,cpzg);
u0 = uexactfn(0,th,phi);

% plotting grid on sphere, based on parameterization
[xp,yp,zp] = sphere(64);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1, yp1, zp1);

% Eplot are cells of matrices which interpolate data onto the plotting grid
n_level = length(a_band);
Eplot = cell(n_level-1, 1);
dx_tmp = dx;
for i = 1:1:n_level-1
    x = (x0:dx_tmp:x1)';
    y = (y0:dx_tmp:y1)';
    z = (z0:dx_tmp:z1)';

    %Eplot{i} = interp3_matrix_band( x, y, z, xp1, yp1, zp1, p, a_band{i} );
    Eplot{i} = interp3_matrix_test( x, y, z, xp1, yp1, zp1, p );
    Eplot{i} = Eplot{i}(:,a_band{i});
    dx_tmp = 2*dx_tmp;
end


%% Time-stepping for the heat equation
Tf = 0.1;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;

gap = 5;
err_plot = zeros(numtimesteps,1);
err_plot1 = zeros(numtimesteps,1);
t_plot = dt:dt:Tf;

u = u0;
figure(1);
disp('starting to do time evolution by matlab')
tic
for kt = 1:numtimesteps
    
    t = kt*dt;
    %tic;
    %[unew flag relres iter] = bicgstab(Mc{1}, u, 1e-10, 20, [], [], u);
    %[unew flag relres iter] = gmres(Mc{1}, u, [], 1e-10, [], [], [], u);
    [unew flag] = gmres(Mc{1}, u, 5, 1e-6);
    u = unew;
    
    %f = u + 0.5*dt*Lc{1}*u;
    %[unew flag relres iter] = gmres(M1, f, [], 1e-6);
    %u = Ec{1}*unew;
    
    %toc;
    
    sphplot = Eplot{1}*u;
    error_sphere_inf = max(abs( uexactfn(t, th_plot, phi_plot) - sphplot )) / max(abs(uexactfn(t, th_plot, phi_plot)));
    err_plot1(kt) = error_sphere_inf;
    %plotting
%     if ( (kt < 5) || (mod(kt,gap) == 0) || (kt == numtimesteps) )
%       [t dt dx error_sphere_inf]
%       % plot value on sphere
%       figure(2); clf;
%       sphplot = reshape(sphplot, size(xp));
%       surf(xp, yp, zp, sphplot);
%       title( ['soln at time ' num2str(t) ', kt= ' num2str(kt)] );
%       xlabel('x'); ylabel('y'); zlabel('z');
%       axis equal; shading interp;
%       colorbar;
%       pause(0.001);
%     end
    
end
time_matlab = toc
disp('done')

close all
u = u0;
V = cell(n_level,1);
F = cell(n_level,1);
figure(1);
disp('starting to do time evolution by multigrid')
tic
for kt = 1:numtimesteps

    t = kt*dt;
    %% implicit backward euler
    % set the initial value of multigrid as the solution of last time step
    
    %tic;
    V{1} = u;
    for i = 2:1:n_level
        V{i} = zeros(size(a_band{i}));
    end
    %F{1} = u + lambda*dt*L1*u;
    F{1} = u;
    %unew = gmg(Mn, Lc, Ec, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w, uexactfn(t,phi_plot), Eplot, MAX);
    unew = helper_vcycle(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w);
    %toc;
    
    u = unew;

    sphplot = Eplot{1}*u;
    error_sphere_inf = max(abs( uexactfn(t, th_plot, phi_plot) - sphplot )) / max(abs(uexactfn(t, th_plot, phi_plot)));
    err_plot(kt) = error_sphere_inf;
    %plotting
%     if ( (kt < 5) || (mod(kt,gap) == 0) || (kt == numtimesteps) )
%       [t dt dx error_sphere_inf]
%       % plot value on sphere
%       figure(2); clf;
%       sphplot = reshape(sphplot, size(xp));
%       surf(xp, yp, zp, sphplot);
%       title( ['soln at time ' num2str(t) ', kt= ' num2str(kt)] );
%       xlabel('x'); ylabel('y'); zlabel('z');
%       axis equal; shading interp;
%       colorbar;
%       pause(0.001);
%     end
    
end
time_mg = toc
disp('done')

figure(3)
plot(t_plot, err_plot1, t_plot,err_plot,'r')
xlabel('time')
ylabel('relative error')
title(['heat equation a sphere  relative error of e^{-2t}cos(\phi),  dx = ', num2str(dx)])
legend('matlab','multigrid')
