%% test geometric multigrid method to solve poisson equation on a circle

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');


x0 = -2;
x1 = 2;

%%
% 2D example on a circle
% Construct a grid in the embedding space

dx = 0.05; % grid size
%dx = 0.05;
dx_coarsest = 0.2;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = x1d_coarsest;
z1d_coarsest = x1d_coarsest;

dy = dx;
dz = dx;

dim = 3;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

n1 = 3;
n2 = 3;

p_f2c = 1;
p_c2f = 1;

w = 1;

cpf = @cpSphere;

has_boundary = false;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg] = ...
    build_mg_cpgrid3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);

n_level = length(a_band);

%disp('building cp matrices ... ')
%[Mc, Lc, Ec] = build_mg_cpmatrix3d(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order);

disp('building cp matrices on each level of V-cycle')
[Mc, Lc, Ec, Eic, Rc, a_iband, a_oband, a_xcp, a_ycp, a_zcp, a_xg, a_yg, a_zg, a_innerInOuter] = ...
    build_mg_ops_and_bands3d(a_band, a_xcp, a_ycp, a_zcp, a_xg, a_yg, a_zg, a_x1d, a_y1d, a_z1d, p, order);
disp('done')

% dangerous in changing lines of E in this way
%for i = 1:1:n_level
%    E = Ec{i};
%	E(a_innerInOuter{i},:) = speye(size(Ec{i},2));
%	%Mc{i} = Lc{i}*E;
%    Ec{i} = E;
%end

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM3d(a_x1d, a_y1d, a_z1d, a_xcp, a_ycp, a_zcp, a_iband, a_bdyg, p_f2c, p_c2f);
disp('done')

disp('extracting diagonal entries of L, making jacobi iteration faster...')
tic;
D = cell(n_level,1);
for i = 1:1:n_level
    [i1,j1,r1] = find(Rc{i});
    Ldiagpad = Rc{i}.*Lc{i};
    D{i} = diag(Ldiagpad(i1,j1));
end
toc;
disp('done') 

%pt = [1,0,0];
%uexactfn = @(th, phi) cos(phi+pi/2);
%rhsfn = @(th, phi, r) -2*cos(phi+pi/2);

%pt = [1,0,0];
%l = 5;
%rhsfn = @(th, phi, r) -2*cos(phi+pi/2) - l*(l+1)*sin(l*th).*sin(phi+pi/2).^l;
%uexactfn = @(th, phi) cos(phi+pi/2) + sin(l*th).*sin(phi+pi/2).^l;

%pt = [1,0,0];
%l = 3;
%rhsfn = @(th, phi, r) -l*(l+1)*sin(l*th).*sin(phi+pi/2).^l;
%uexactfn1 = @(th, phi) sin(l*th).*sin(phi+pi/2).^l;
%uexactfn = @(th, phi) uexactfn1(th, phi) - uexactfn1(0,0);

pt = [0,0,1];
rhsfn = @(th, phi, r) ( -30*cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1) );
uexactfn1 = @(th, phi) cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1);
uexactfn = @(th, phi) uexactfn1(th, phi) - uexactfn1(0,pi/2);


%% building E_plot for purpose of plotting and debug
% plotting grid on circle, using theta as a parameterization
cpxg = a_xcp{1};
cpyg = a_ycp{1};
cpzg = a_zcp{1};
[th, phi, r] = cart2sph(cpxg,cpyg,cpzg);

u0 = uexactfn(th, phi);

[xp, yp, zp] = sphere(64);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);

n_level = length(a_band);
Eplot = cell(n_level,1);
dx_tmp = dx;
for i = 1:1:n_level
    x = (x0:dx_tmp:x1)';
    y = x;
    z = x;
    
    Eplot{i} = interp3_matrix_test( x, y, z, xp1, yp1, zp1, p );
    Eplot{i} = Eplot{i}(:,a_iband{i});
    dx_tmp = 2*dx_tmp;
end

%% Setting up right hand side

uexact = uexactfn(th_plot, phi_plot);

uexact_debug = cell(n_level,1);
for i = 1:1:n_level
    [th, phi, r] = cart2sph(a_xcp{i}, a_ycp{i}, a_zcp{i});
    uexact_debug{i} = uexactfn(th, phi);
end

disp('building right hand side and allocate space for solution ... ')
[V, F] = helper_set_rhs3d(a_xcp, a_ycp, a_zcp, rhsfn, 1);
disp('done')

disp('making the problem to solve poisson equation on a circle well-posedness ... ')

% pt is a point whose value we want to specify as the same at each level of
% V-Cycle.

% Following line of code seems a bit confusing because a circle does not have a
% boundary, and the logical variable 'has_boundary' is set to be 'false',
% but in order to make the problem wellposedness, we do something with the
% coefficient matrix just as when dealing with Neumann Boundary Conditions,
% so we add the following line:
% [Mc, Lc, Ec] = app_bnd3d(Mc, Lc, Ec, F, a_xcp, a_ycp, a_zcp, a_bdyg, 'neumann');
 [Mc, Lc, Ec] = app_bnd3d_2b(Mc, Lc, Ec, F, a_innerInOuter, a_xg, a_yg, a_zg, a_bdyg, pt, 'neumann');
disp('done')


% circplot = cell(n_level-1,1);
% error_inf_matlab = zeros(n_level-1,1);
% res_matlab = zeros(n_level,1);
% for i = 1:1:n_level-1
% unew = Mc{i} \ F{i};
% 
% % By increasing the maximal number of iterations, bicgstab will converge;
% % however as grid become finer, number of iterations will increase; and for
% % dx less or equal to 0.00625, bicgstab does not converge after 2000
% % iterations
% % [unew flag] = bicgstab(Mc{i}, F{i}, 1e-10, 200);
% 
% % gmres seems not converge for this problem
% % [unew flag] = gmres(Mc{i}, F{i}, 10, 1e-10, 200);
% 
% circplot{i} = Eplot{i}*unew;
% error_inf_matlab(i) = max(abs( uexact - circplot{i} ));
% res_matlab(i) = norm(Eplot{i}*(F{i} - Mc{i}*unew),inf);
% 
% end

i = n_level;
tic;
unew = Mc{i} \ F{i};
t_matlab1 = toc
circplot = Eplot{i}*unew;
err_inf_matlab1 = max(abs(uexact-circplot));

i = n_level-1;
tic;
unew = Mc{i} \ F{i};
t_matlab2 = toc
circplot = Eplot{i}*unew;
err_inf_matlab2 = max(abs(uexact-circplot));


MAX = 50;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
for start = 1:1:n_level-1
%     [thg rg] = cart2pol(a_xcp{start},a_ycp{start});
%     V{start} = v_initial_fn(thg);
%    F{start} = zeros(size(F{start}));
%    for i = start:1:n_level
%        V{i} = zeros(size(V{i}));
%    end
    tic;
    %[umg err_inf(start,:) res(start,:)] = gmg_M(Mc, Lc, Eic, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, start, w, uexact, Eplot, MAX); 
    [umg err_inf(start,:) res(start,:)] = gmg_2b(Mc, Lc, Ec, Eic, D, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, start, w, uexact, Eplot, MAX);
    t_mg = toc
end

sphplot = Eplot{n_level-1}*umg;
%sphplot = Eplot{1}*F{1};
%sphplot = uexact;
figure(1)
sphplot = reshape(sphplot, size(xp));
surf(xp, yp, zp, sphplot);
title( 'numerical solution of poisson equation on a sphere with exact solution as cos(\phi) ');
xlabel('x'); ylabel('y'); zlabel('z');
axis equal; shading interp;
colorbar;

figure(2)
% rep_res_matlab = repmat(res_matlab, 1, 2);
% xx = [0 7];
% semilogy(xx,rep_res_matlab(1,:),'b',xx,rep_res_matlab(2,:),'r',xx,rep_res_matlab(3,:),'c', ...
%          xx,rep_res_matlab(4,:),'k');
% hold on

n = 1:MAX;
%semilogy(n,res(1,:),'.-', n,res(2,:),'r.-', n,res(3,:),'c^-', n,res(4,:),'k-s');
%legend('N=80','N=40','N=20','N=10')
%semilogy(n,res(1,:),'.-', n,res(2,:),'r.-', n,res(3,:),'c^-');
%legend('N=40','N=20','N=10')
semilogy(n,res(1,:),'.-', n,res(2,:),'r.-');
legend('N=20','N=10')
xlabel('number of vcyles')
ylabel('|residual|_{\infty}')
title(['cos(\phi) with p=', num2str(p), ',  res = E*(f-L*v)'])

% plot error of matlab and error of different number of vcycles
figure(3)

%err_inf_matlab = cell2mat(error_inf_matlab);
%rep_err_inf_matlab = repmat(error_inf_matlab,1,2);
%xx = [0 7];
%semilogy(xx,rep_err_inf_matlab(1,:),'b',xx,rep_err_inf_matlab(2,:),'r',xx,rep_err_inf_matlab(3,:),'c', ...
%         xx,rep_err_inf_matlab(4,:),'k');
%semilogy(xx,rep_err_inf_matlab(1,:),'b',xx,rep_err_inf_matlab(2,:),'r');
%hold on

n = 1:MAX;
%semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-', n,err_inf(4,:), 'ko-');
%legend('N=80','N=40','N=20','N=10')
%semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-');
%legend('N=40','N=20','N=10')
semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
legend('N=20','N=10')
xlabel('number of vcyles')
ylabel('|error|_{\infty}')
title(['cos(\phi) with p=', num2str(p)])
