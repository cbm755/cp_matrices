%% test geometric multigrid method to solve poisson equation on a sphere

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

dx = 0.1; % grid size
%dx = 0.025;
%dx = 0.00625;  % lots of memory...
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

cpf = @cpSphere;

has_boundary = false;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg] = ...
    build_mg_cpgrid3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);

n_level = length(a_band);

shift = 1;
C = 1.5;
m = 1;

% syms theta phi r
% beta_fun = cos(phi+pi/2)+C;
% uexactfn = cos(phi+pi/2);
% rhsfn = 1/sin(phi+pi/2)*( diff( beta_fun/sin(phi+pi/2)*diff(uexactfn,theta), theta) + ...
%                      diff( beta_fun*sin(phi+pi/2)*diff(uexactfn,phi), phi) );
%                 
% beta_fun = matlabFunction(beta_fun,'vars',{theta,phi,r}); 
% uexactfn = matlabFunction(uexactfn,'vars',{theta,phi,r});
% rhsfn = matlabFunction(rhsfn,'vars',{theta,phi,r}); 
% rhsfn = @(theta, phi, r) rhsfn(theta,phi,r) - shift*uexactfn(theta,phi,r);

beta_fun = @(theta,phi,r) cos(m*(phi+pi/2))+C;
uexactfn = @(theta,phi,r) cos(phi+pi/2);
rhsfn = @(theta,phi,r) -2*cos(phi+pi/2).*beta_fun(theta,phi,r) + m*sin(phi+pi/2).*sin(m*(phi+pi/2)) - shift*uexactfn(theta,phi,r);

% beta_fun = @(theta,phi,r) ones(size(theta));
% uexactfn = @(theta,phi,r) cos(phi+pi/2);
% rhsfn = @(theta,phi,r) -2*cos(phi+pi/2) -shift*uexactfn(theta,phi,r);
                 

disp('building cp matrices ... ')
%[Mc, Lc, Ec] = build_mg_cpmatrix3d(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order);
Mc = cell(n_level,1);
Lc = cell(n_level,1);
Ec = cell(n_level,1);

for i = 1:1:n_level
    E1 = interp3_matrix(a_x1d{i},a_y1d{i},a_z1d{i},a_xcp{i},a_ycp{i},a_zcp{i},1);
    E1 = E1(:,a_band{i});
    
    Ec{i} = interp3_matrix(a_x1d{i},a_y1d{i},a_z1d{i},a_xcp{i},a_ycp{i},a_zcp{i},p);
    Ec{i} = Ec{i}(:,a_band{i});

    % If the coefficients were constant, just call 'laplacian_2d_matrix' then done;
	% for the variable coefficient case, we need a few more work to do...
    
    % First compute the forward and backward difference in the x and y
    % direction
    [Dxb, Dxf, Dyb, Dyf, Dzb, Dzf] = firstderiv_upw1_3d_matrices(a_x1d{i},a_y1d{i},a_z1d{i},a_band{i});
    
    % Then compute the coefficients: beta_{i+1/2,j}, beta_{i-1/2,j},
    % beta_{i,j+1/2}, and beta_{i,j-1/2}.
    [theta, phi, r] = cart2sph(a_xcp{i}, a_ycp{i}, a_zcp{i});
    beta = beta_fun(theta, phi, r);
    [Axb, Axf, Ayb, Ayf, Azb, Azf] = average_upw1_3d_matrices(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_band{i});
    coef_xb = Axb*beta; 
    coef_xf = Axf*beta;
    coef_yb = Ayb*beta;
    coef_yf = Ayf*beta;
    coef_zb = Azb*beta;
    coef_zf = Azf*beta;
    % We just want to multiply the value of beta with each corresponding row of 
    % the difference matrices, but making the vector beta into a sparse
    % diagonal matrix, and then multiply the difference matrices might be
    % faster; another choice would be 'bsxfun'; better way to do this?
	s = length(a_band{i});
    diag_coef_xb = spdiags(coef_xb,0,s,s); 
    diag_coef_xf = spdiags(coef_xf,0,s,s);
    diag_coef_yb = spdiags(coef_yb,0,s,s);
    diag_coef_yf = spdiags(coef_yf,0,s,s);
    diag_coef_zb = spdiags(coef_zb,0,s,s);
    diag_coef_zf = spdiags(coef_zf,0,s,s);
	
	Lc{i} = (diag_coef_xf*Dxf - diag_coef_xb*Dxb)/a_dx{i} + ...
	        (diag_coef_yf*Dyf - diag_coef_yb*Dyb)/a_dx{i} + ...
            (diag_coef_zf*Dzf - diag_coef_zb*Dzb)/a_dx{i};

    I = speye(size(Lc{i}));
    
    % E3*L is better than E1*L
    Mc{i} = E1*Lc{i}- 2*dim/a_dx{i}^2*(I-Ec{i});
	
    % Because of the shift, E3*L*E3 also works
    %Mc{i} = Ec{i}*Lc{i}*Ec{i};
    
    Lc{i} = Lc{i} - shift*I;
    Mc{i} = Mc{i} - shift*I;
end

disp('done')

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM3d(a_x1d, a_y1d, a_z1d, a_xcp, a_ycp, a_zcp, a_band, a_bdyg, p_f2c, p_c2f);
disp('done')
 

%% building E_plot for purpose of plotting and debug
% plotting grid on circle, using theta as a parameterization
cpxg = a_xcp{1};
cpyg = a_ycp{1};
cpzg = a_zcp{1};
[th, phi, r] = cart2sph(cpxg,cpyg,cpzg);

u0 = uexactfn(th, phi);

[xp, yp, zp] = sphere(256);
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
%[V, F] = helper_set_rhs3d(a_xg, a_yg, a_zg, rhsfn, 1);
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
%unew = agmg(Mc{i},F{i},[],1e-10,1000);
t_matlab1 = toc
circplot = Eplot{i}*unew;
err_inf_matlab1 = max(abs(uexact-circplot));

i = n_level-1;
tic;
unew = Mc{i} \ F{i};
%unew = agmg(Mc{i},F{i},[],1e-10,1000);
t_matlab2 = toc
circplot = Eplot{i}*unew;
err_inf_matlab2 = max(abs(uexact-circplot));

i = n_level-2;
tic;
unew = Mc{i} \ F{i};
%unew = agmg(Mc{i},F{i},[],1e-10,1000);
t_matlab3 = toc
circplot = Eplot{i}*unew;
err_inf_matlab3 = max(abs(uexact-circplot));
% 
% i = n_level-3;
% tic;
% unew = Mc{i} \ F{i};
% %unew = agmg(Mc{i},F{i},[],1e-10,1000);
% t_matlab4 = toc
% circplot = Eplot{i}*unew;
% err_inf_matlab4 = max(abs(uexact-circplot));

% i = n_level-4;
% tic;
% unew = Mc{i} \ F{i};
% %unew = agmg(Mc{i},F{i},[],1e-10,1000);
% t_matlab5 = toc
% circplot = Eplot{i}*unew;
% err_inf_matlab5 = max(abs(uexact-circplot));

% i = n_level-5;
% tic;
% unew = Mc{i} \ F{i};
% %unew = agmg(Mc{i},F{i},[],1e-10,1000);
% t_matlab6 = toc
% circplot = Eplot{i}*unew;
% err_inf_matlab6 = max(abs(uexact-circplot));

MAX = 50;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
R = [];
umg = cell(n_level,1);
for start = 1:1:n_level-1
%     [thg rg] = cart2pol(a_xcp{start},a_ycp{start});
%     V{start} = v_initial_fn(thg);
%    F{start} = zeros(size(F{start}));
%    for i = start:1:n_level
%        V{i} = zeros(size(V{i}));
%    end
    tic;
    [umg{start}, err_inf(start,:), res(start,:)] = gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, Eplot, MAX);
    t_mg = toc
end

err_inf = err_inf(end:-1:1,:);
res = res(end:-1:1,:);

i = n_level-2;
sphplot = Eplot{i}*umg{i};
%sphplot = Eplot{1}*F{1};
%sphplot = uexact;

error_inf_matlab = [0.098803124412236;
   0.022296792714046;
   0.005656587868815;
   0.001389990501488;
   0.000349262823944
   nan];

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
elseif n_level == 6
semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', n,err_inf(4,:), 'k-s', ...
         n,err_inf(5,:),'c^-');
legend('N=5','N=10','N=20','N=40','N=80')
elseif n_level == 5
semilogy(n,res(1,:),'o--', n,res(2,:),'r*--', n,res(3,:),'g+--', n,res(4,:),'k-s');
legend('N=5','N=10','N=20','N=40')
elseif n_level == 4
semilogy(n,res(1,:),'o--', n,res(2,:),'r*--', n,res(3,:),'g+--');
legend('N=5','N=10','N=20')
end
% semilogy(n,res(1,:),'.-', n,res(2,:),'r.-');
% legend('N=20','N=10')
xlabel('\fontsize{15} number of vcyles')
ylabel('\fontsize{15} ||f^h-A^hu^h||_{\infty}/||f^h||_{\infty}')
%title('\fontsize{15} relative residuals')
set(gca,'Fontsize',12)
xlim([0 10])

% plot error of matlab and error of different number of vcycles
figure(3)
n = 1:MAX;
n = n-1;
if n_level == 7
semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', n,err_inf(4,:), 'k-s', ...
         n,err_inf(5,:),'c^-',n,err_inf(6,:),'m-d');
legend('N=5','N=10','N=20','N=40','N=80','N=160')
elseif n_level == 6
semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', n,err_inf(4,:), 'k-s', ...
         n,err_inf(5,:),'c^-');
legend('N=5','N=10','N=20','N=40','N=80')
elseif n_level == 5
semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', n,err_inf(4,:), 'k-s');
legend('N=5','N=10','N=20','N=40')
elseif n_level == 4
semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--');
legend('N=5','N=10','N=20')
end

hold on,
%rep_err_inf_matlab = repmat(error_inf_matlab,1,2);
%xx = [0 7];
%semilogy(xx,rep_err_inf_matlab(1,:),'--',xx,rep_err_inf_matlab(2,:),'r--',xx,rep_err_inf_matlab(3,:),'g', ...
%         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'c',xx,rep_err_inf_matlab(6,:),'m');

% semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
% legend('N=20','N=10')
xlabel('\fontsize{15}number of vcyles')
% ylabel('\fontsize{15}|error|_{\infty}')
ylabel('\fontsize{15} ||u^h-u||_{\infty}/||u||_{\infty}')
%title('\fontsize{15} relative errors')
%title('\fontsize{15} u=cos(3\theta)sin(\phi)^3(9cos(\phi)^2-1)')
%title('\fontsize{15}cos(phi), Dirichlet B.C.s')
set(gca,'Fontsize',12)
xlim([0 10])
