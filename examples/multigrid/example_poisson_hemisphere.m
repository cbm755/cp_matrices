%% test geometric multigrid method to solve poisson equation on a hemisphere

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');

% add notay amg
addpath('/scratch/cheny1/opt/AGMG_3.1.1/Matlab')

x0 = -3;
x1 = 3;

%%
% 2D example on a circle
% Construct a grid in the embedding space

dx = 0.00625; % grid size
%dx = 0.05;
dx_coarsest = 0.4;   % coarsest grid size
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

cpf1 = @cpHemisphere;  paramf = @paramHemisphere; cpf = @(x,y,z) cpbar_3d(x,y,z,cpf1);

has_boundary = true;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg] = ...
    build_mg_cpgrid3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);
disp('done')

n_level = length(a_band);

disp('building cp matrices ... ')
Mc = cell(n_level,1);
Lc = cell(n_level,1);
Ec = cell(n_level,1);
for i = 1:1:n_level
   ddx = a_x1d{i}(2) - a_x1d{i}(1);
   Ec{i} = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, p);
   Ec{i} = Ec{i}(:, a_band{i});
   Lc{i} = laplacian_3d_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i});
   E = interp3_matrix(a_x1d{i}, a_y1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, 1);
   E = E(:,a_band{i});
%      bdy = logical(a_bdyg{i});
%      E(bdy,:) = - E(bdy,:);
%      Ec{i}(bdy,:) = - Ec{i}(bdy,:);
   Mc{i} = E*Lc{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});
   %Mc{i} = Lc{i}*Ec{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});
end
shift = 1;
for i = 1:1:n_level
    Mc{i} = Mc{i} - shift*speye(size(Mc{i}));
    Lc{i} = Lc{i} - shift*speye(size(Lc{i}));
end

%[Mc, Lc, Ec] = build_mg_cpmatrix3d(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order);
disp('done')

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM3d(a_x1d, a_y1d, a_z1d, a_xcp, a_ycp, a_zcp, a_band, a_bdyg, p_f2c, p_c2f);
disp('done')
 
R = cell(n_level,1);
for i = 1:1:n_level
    R{i} = a_xg{i}.^2 + a_yg{i}.^2 + a_zg{i}.^2;
end

% For Dirichlet boundary condition

%   uexactfn = @(th, phi) cos(phi+pi/2);
%   rhsfn = @(th, phi, r) -2*cos(phi+pi/2)./(r.^2);

% For Nuemann boundary condition

% pt is a point whose value we want to specify as the same at each level of
% V-Cycle.
%pt = [1,0,0];


% shift version
uexactfn = @(th, phi,r) sin(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1);
rhsfn = @(th,phi,r) -shift*uexactfn(th,phi) + ( -30*sin(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1) );

% converge fast
%rhsfn = @(th, phi, r) ( -30*sin(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1) );
%uexactfn1 = @(th, phi) sin(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1);
%uexactfn = @(th, phi) uexactfn1(th, phi) - uexactfn1(0,0);

% converge fast
%l = 5;
%rhsfn = @(th, phi, r) ( -l*(l+1)*sin(l*th).*sin(phi+pi/2).^l ) ./ (r.^2);
%uexactfn1 = @(th, phi) sin(l*th).*sin(phi+pi/2).^l;
%uexactfn = @(th, phi) uexactfn1(th, phi) - uexactfn1(0,0);


% converge not that fast and need to specify another point:
%  pt = [1 0 0];
 
%  uexactfn = @(th, phi) cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1);
%  rhsfn = @(th, phi, r) ( -30*cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1) ) - shift*uexactfn(th,phi);
  
%  uexactfn1 = @(th, phi) cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1);
%  uexactfn = @(th, phi) uexactfn1(th, phi) - uexactfn1(0,0);

%  pt = [0 0 1];
% 
%  rhsfn = @(th, phi, r) ( -30*cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1) );
%  uexactfn1 = @(th, phi) cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1);
%  uexactfn = @(th, phi) uexactfn1(th, phi) - uexactfn1(0,pi/2) + rhsfn(0,pi/2); 


% converge not that fast
%l = 5;
%rhsfn = @(th, phi, r) -l*(l+1)*cos(l*th).*sin(phi+pi/2).^l;
%uexactfn = @(th, phi) cos(l*th).*sin(phi+pi/2).^l;

%uexactfn1 = @(th, phi) cos(l*th).*sin(phi+pi/2).^l;
%uexactfn = @(th, phi) uexactfn1(th, phi) - uexactfn1(0,0);


%% building E_plot for purpose of plotting and debug
% plotting grid on circle, using theta as a parameterization
cpxg = a_xcp{1};
cpyg = a_ycp{1};
cpzg = a_zcp{1};
[th, phi, r] = cart2sph(cpxg,cpyg,cpzg);

u0 = uexactfn(th, phi);

[xp, yp, zp] = paramHemisphere(256);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);

n_level = length(a_band);
Eplot = cell(n_level,1);
dx_tmp = dx;
for i = 1:1:n_level
    x = (x0:dx_tmp:x1)';
    y = x;
    z = x;
    
    Eplot{i} = interp3_matrix( x, y, z, xp1, yp1, zp1, p );
    Eplot{i} = Eplot{i}(:,a_band{i});
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
% for i = 1:1:n_level
%     bdy = logical(a_bdyg{i});
%     F{i}(bdy,:) = - F{i}(bdy,:);
% end
%[V, F] = helper_set_rhs3d(a_xg, a_yg, a_zg, rhsfn, 1);
disp('done')

disp('deal with boundary conditions ... ')

% tol = 1e-10;
 Mc1 = Mc;
 Lc1 = Lc;
 F1 = F;
% i = n_level;
% %pt1 = [-0.2 1.4 0.2];
% pt1 = [0 0 1];
% j = abs(a_xg{i}-pt1(1)) < tol & abs(a_yg{i}-pt1(2)) < tol & abs(a_zg{i}-pt1(3)) < tol;
% Lc1{i}(j,:) = sum(Ec{i});
% Mc1{i}(j,:) = sum(Ec{i});
% F1{i}(j) = 0;

pt1 = [0.2 0.2 1];

% [Mc, Lc, Ec, F] = app_bnd3d(Mc, Lc, Ec, F, a_xg, a_yg, a_zg, a_bdyg, pt1, 'dirichlet');
% [Mc, Lc, Ec, F] = app_bnd3d(Mc, Lc, Ec, F, a_xg, a_yg, a_zg, a_bdyg, pt1, 'neumann');
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

TOL = 1e-10;
MAXIT = 1000;
VERBOSE = 0;

i = n_level;
tic;
unew = Mc{i} \ F{i};
%unew = agmg(Mc{i},F{i},[],1e-10,1000);
t_matlab0 = toc
circplot = Eplot{i}*unew;
err_inf_matlab0 = max(abs(uexact-circplot)) / norm(uexact,inf);

i = n_level-1;
tic;
unew = Mc{i} \ F{i};
%unew = agmg(Mc{i},F{i},[],1e-10,1000);
t_matlab1 = toc
circplot = Eplot{i}*unew;
err_inf_matlab1 = max(abs(uexact-circplot)) / norm(uexact,inf);

i = n_level-2;
tic;
unew = Mc{i} \ F{i};
%unew = agmg(Mc{i},F{i},[],1e-10,1000);
t_matlab2 = toc
circplot = Eplot{i}*unew;
err_inf_matlab2 = max(abs(uexact-circplot)) / norm(uexact,inf);

% i = n_level-3;
% tic;
% unew = Mc{i} \ F{i};
% %unew = agmg(Mc{i},F{i},[],1e-10,1000);
% t_matlab3 = toc
% circplot = Eplot{i}*unew;
% err_inf_matlab3 = max(abs(uexact-circplot)) / norm(uexact,inf);

%% t = 11
%i = n_level-4;
%tic;
%unew = Mc{i} \ F{i};
%%unew = agmg(Mc{i},F{i},[],1e-10,1000);
%t_matlab4 = toc
%circplot= Eplot{i}*unew;
%err_inf_matlab4 = max(abs(uexact-circplot)) / norm(uexact,inf);
%% t = 77
% i = n_level-5;
% tic;
% unew = Mc{i} \ F{i};
% %unew = agmg(Mc{i},F{i},[],TOL,MAXIT);
% t_matlab5 = toc
% circplot= Eplot{i}*unew;
% err_inf_matlab5 = max(abs(uexact-circplot)) / norm(uexact,inf);
%% t = ... 
% i = n_level-6;
% tic;
% unew = Mc{i} \ F{i};
% t_matlab6 = toc
% circplot= Eplot{i}*unew;
% err_inf_matlab6 = max(abs(uexact-circplot)) / norm(uexact,inf);

%% relative errors of matlab's backslash, Dirichelt B.C.s:
% error_inf_matlab = [ 0.006190353851154;
%     0.001156657183076;
%     2.455153673996957e-04;
%     5.582573195050955e-05;
%     1.320532238113792e-05;
%     3.214392619188544e-06];


%% absolute errors of matlab's backslash, Neumann B.C.s:
error_inf_matlab = [ 0.0988094245253703;
        0.0223065773809175;
         0.005658536036847;
       0.00139026894892236;
      0.000349291551987055;
      8.71744968649444e-05];


MAX = 50;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
umg = cell(n_level,1);
for start = 1:1:n_level-1
%     [thg rg] = cart2pol(a_xcp{start},a_ycp{start});
%     V{start} = v_initial_fn(thg);
%    F{start} = zeros(size(F{start}));
%    for i = start:1:n_level
%        V{i} = zeros(size(V{i}));
%    end
    tic;
    %[umg err_inf(start,:) res(start,:)] = gmg(Mc1, Lc1, Ec, V, F1, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, Eplot, MAX);
    [umg{start} err_inf(start,:) res(start,:)] = gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, Eplot, MAX);
    t_mg = toc
end
err_inf = err_inf(end:-1:1,:);
res = res(end:-1:1,:);

i = n_level-2;
sphplot = Eplot{i}*umg{i};
%sphplot = Eplot{1}*F{1};
%sphplot = uexact;

figure(1)
sphplot = reshape(sphplot, size(xp));
surf(xp, yp, zp, sphplot);
title( '\fontsize{15} cos(3\theta)sin(\phi)^3(9cos(\phi)^2-1) ');
xlabel('\fontsize{15}x'); ylabel('\fontsize{15}y'); zlabel('\fontsize{15}z');
axis equal; shading interp;
colorbar;
set(gca,'Fontsize',12)

figure(2)
% rep_res_matlab = repmat(res_matlab, 1, 2);
% xx = [0 7];
% semilogy(xx,rep_res_matlab(1,:),'b',xx,rep_res_matlab(2,:),'r',xx,rep_res_matlab(3,:),'c', ...
%          xx,rep_res_matlab(4,:),'k');
% hold on

n = 1:MAX;
n = n-1;
if n_level == 7
semilogy(n,res(1,:),'o--', n,res(2,:),'r*--', n,res(3,:),'g+--', n,res(4,:),'k-s', ...
         n,res(5,:),'c^-',n,res(6,:),'m-d');
legend('N=5','N=10','N=20','N=40','N=80','N=160')
elseif n_level == 5
semilogy(n,res(1,:),'o--', n,res(2,:),'r*--', n,res(3,:),'g+--', n,res(4,:),'k-s');
legend('N=5','N=10','N=20','N=40')
elseif n_level == 4
semilogy(n,res(1,:),'o--', n,res(2,:),'r*--', n,res(3,:),'g+--');
legend('N=5','N=10','N=20')
end
% semilogy(n,res(1,:),'.-', n,res(2,:),'r.-');
% legend('N=20','N=10')
xlabel('number of vcyles')
ylabel('|residual|_{\infty}')
title(['cos(\phi) with p=', num2str(p), ',  res = E*(f-L*v)'])

% plot error of matlab and error of different number of vcycles
figure(3)
n = 1:MAX;
n = n-1;
if n_level == 7
semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', n,err_inf(4,:), 'k-s', ...
         n,err_inf(5,:),'c^-',n,err_inf(6,:),'m-d');
legend('N=5','N=10','N=20','N=40','N=80','N=160')
elseif n_level == 5
semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', n,err_inf(4,:), 'k-s');
legend('N=5','N=10','N=20','N=40')
elseif n_level == 4
semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--');
legend('N=5','N=10','N=20')
end

hold on,
rep_err_inf_matlab = repmat(error_inf_matlab,1,2);
xx = [0 7];
semilogy(xx,rep_err_inf_matlab(1,:),'--',xx,rep_err_inf_matlab(2,:),'r--',xx,rep_err_inf_matlab(3,:),'g', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'c',xx,rep_err_inf_matlab(6,:),'m');

% semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
% legend('N=20','N=10')
xlabel('\fontsize{15}number of vcyles')
% ylabel('\fontsize{15}|error|_{\infty}')
ylabel('\fontsize{15} ||u^h-u||_{\infty}/||u||_{\infty}')
title('\fontsize{15}cos(3\theta)sin(\phi)^3(9cos(\phi)^2-1), Neumann B.C.s')
%title('\fontsize{15}cos(phi), Dirichlet B.C.s')
set(gca,'Fontsize',12)
