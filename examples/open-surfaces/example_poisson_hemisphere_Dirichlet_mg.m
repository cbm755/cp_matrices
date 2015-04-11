%% Geometric Multigrid Poisson equation on a hemisphere with Dirichlet B.C.

%% Using cp_matrices
x0 = -4;
x1 = 4;

%%
% Construct a grid in the embedding space

%dx = 0.0125 % minimum
dx = 0.0125; % grid size
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

n1 = 5;
n2 = 5;

p_f2c = 1;
p_c2f = 1;

w = 1;

R = sqrt(2);
cpf = @(x,y,z) cpHemisphere(x,y,z,R);  paramf = @paramHemisphere;

has_boundary = true;

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
   Ec{i} = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, p);
   Ec{i} = Ec{i}(:, a_band{i});
   Lc{i} = laplacian_3d_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i});
   E = interp3_matrix(a_x1d{i}, a_y1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, 1);
   E = E(:,a_band{i});
   Mc{i} = E*Lc{i} - 2*dim/a_dx{i}^2*(speye(size(E))-Ec{i});
end
shift = 0;
for i = 1:1:n_level
    Mc{i} = Mc{i} - shift*speye(size(Mc{i}));
    Lc{i} = Lc{i} - shift*speye(size(Lc{i}));
end

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM3d(a_x1d, a_y1d, a_z1d, a_xcp, a_ycp, a_zcp, a_band, a_bdyg, p_f2c, p_c2f);
disp('done')

% For Dirichlet boundary condition

%   uexactfn = @(th, phi) cos(phi+pi/2);
%   rhsfn = @(th, phi, r) -2*cos(phi+pi/2)./(r.^2);

% shift version
uexactfn = @(th, phi) cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1);
rhsfn = @(th, phi, r) ( -30*cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1) )./r.^2 - shift*uexactfn(th,phi);

% converge fast
%rhsfn = @(th, phi, r) ( -30*sin(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1) );
%uexactfn1 = @(th, phi) sin(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1);
%uexactfn = @(th, phi) uexactfn1(th, phi) - uexactfn1(0,0);

% converge fast
%l = 5;
%rhsfn = @(th, phi, r) ( -l*(l+1)*sin(l*th).*sin(phi+pi/2).^l ) ./ (r.^2);
%uexactfn1 = @(th, phi) sin(l*th).*sin(phi+pi/2).^l;
%uexactfn = @(th, phi) uexactfn1(th, phi) - uexactfn1(0,0);

%% building E_plot for purpose of plotting and debug
% plotting grid on circle, using theta as a parameterization

[xp, yp, zp] = paramHemisphere(256,R);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);

n_level = length(a_band);
Eplot = cell(n_level,1);
dx_tmp = dx;
for i = 1:1:n_level
    Eplot{i} = interp3_matrix( a_x1d{i}, a_y1d{i}, a_z1d{i}, xp1, yp1, zp1, p, a_band{i} );
end

%% Setting up right hand side
uexact = uexactfn(th_plot, phi_plot);

uexact_debug = cell(n_level,1);
for i = 1:1:n_level
    [th, phi, r] = cart2sph(a_xcp{i}, a_ycp{i}, a_zcp{i});
    uexact_debug{i} = uexactfn(th, phi);
end

disp('setting up rhs and allocate spaces for solns')
F = cell(n_level,1);
V = cell(n_level,1);
for i = 1:1:n_level
    % Caution: use cart2sph instead of cart2pol!!
    [th, phi, r] = cart2sph(a_xcp{i},a_ycp{i},a_zcp{i});
    F{i} = rhsfn(th,phi,r);
    bdyg = logical(a_bdyg{i});
    F{i}(bdyg) = uexactfn(th(bdyg), phi(bdyg));
    V{i} = zeros(size(F{i}));
end

disp('buidling matrices to deal with boundary conditions ... ')
E_out_out = cell(n_level,1);
E_out_in = cell(n_level,1); 
a_Ebar = cell(n_level,1);
a_Edouble = cell(n_level,1);
a_Etriple = cell(n_level,1);
for i = 1:1:n_level
    x1d = a_x1d{i}; y1d = a_y1d{i}; z1d = a_z1d{i}; band = a_band{i};
    I = speye(size(Lc{i}));
    bdy = logical(a_bdyg{i});
    xg_bar = 2*a_xcp{i}(bdy) - a_xg{i}(bdy);
    yg_bar = 2*a_ycp{i}(bdy) - a_yg{i}(bdy);
    zg_bar = 2*a_zcp{i}(bdy) - a_zg{i}(bdy);
    [cpx_bar,cpy_bar,cpz_bar] = cpf(xg_bar,yg_bar,zg_bar);
    Ebar = interp3_matrix(x1d,y1d,z1d,cpx_bar,cpy_bar,cpz_bar,p,band);
    xg_double = 2*xg_bar - a_xcp{i}(bdy);
    yg_double = 2*yg_bar - a_ycp{i}(bdy);
    zg_double = 2*zg_bar - a_zcp{i}(bdy); 
    [cpx_double, cpy_double, cpz_double] = cpf(xg_double,yg_double,zg_double);
    Edouble = interp3_matrix(x1d,y1d,z1d,cpx_double,cpy_double,cpz_double,p,band);
    xg_triple = 2*xg_double - xg_bar;
    yg_triple = 2*yg_double - yg_bar;
    zg_triple = 2*zg_double - zg_bar;
    [cpx_triple, cpy_triple, cpz_triple] = cpf(xg_triple,yg_triple,zg_triple);
    Etriple = interp3_matrix(x1d,y1d,z1d,cpx_triple,cpy_triple,cpz_triple,p,band);
    M_bdy = (I(bdy,:) + Ebar)/2;
    %M_bdy = (I(bdy,:) + 3*Ebar - Edouble) / 3;
    %M_bdy = (I(bdy,:) + 6*Ebar - 4*Edouble + Etriple) / 4;
    E_out_out{i} = M_bdy(:,bdy);
    E_out_in{i} = M_bdy(:,~bdy);
    Mc{i}(bdy,:) = M_bdy; 
    a_Ebar{i} = Ebar;
    a_Edouble{i} = Edouble;
    a_Etriple{i} = Etriple;
end 

TOL = 1e-10;
MAXIT = 1000;
VERBOSE = 0;

i = n_level;
tic;
unew = Mc{i} \ F{i};
%unew = agmg(Mc{i},F{i},[],1e-10,1000);
t_matlab0 = toc
sphplot = Eplot{i}*unew;
err_inf_matlab0 = max(abs(uexact-sphplot)) / norm(uexact,inf);

i = n_level-1;
tic;
unew = Mc{i} \ F{i};
%unew = agmg(Mc{i},F{i},[],1e-10,1000);
t_matlab1 = toc
sphplot = Eplot{i}*unew;
err_inf_matlab1 = max(abs(uexact-sphplot)) / norm(uexact,inf);

i = n_level-2;
tic;
unew = Mc{i} \ F{i};
%unew = agmg(Mc{i},F{i},[],1e-10,1000);
t_matlab2 = toc
circplot = Eplot{i}*unew;
err_inf_matlab2 = max(abs(uexact-circplot)) / norm(uexact,inf);

i = n_level-3;
tic;
unew = Mc{i} \ F{i};
%unew = agmg(Mc{i},F{i},[],1e-10,1000);
t_matlab3 = toc
circplot = Eplot{i}*unew;
err_inf_matlab3 = max(abs(uexact-circplot)) / norm(uexact,inf);

error_inf_matlab = [err_inf_matlab0; err_inf_matlab1; err_inf_matlab2; err_inf_matlab3]
order_matlab = log(error_inf_matlab(1:end-1)./error_inf_matlab(2:end))/log(2)


%% t = 11
% i = n_level-4;
% tic;
% unew = Mc{i} \ F{i};
% %unew = agmg(Mc{i},F{i},[],1e-10,1000);
% t_matlab4 = toc
% circplot= Eplot{i}*unew;
% err_inf_matlab4 = max(abs(uexact-circplot)) / norm(uexact,inf);
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

MAX = 50;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
umg = cell(n_level,1);
for start = 1:1:n_level-1
    V{start} = zeros(size(F{start}));
    %V{start} = ones(size(F{start}));
    %V{start} = rand(size(F{start})) - 0.5;
    for i = start+1:1:n_level
        V{i} = zeros(size(F{i}));
    end
    [umg{start}, err_inf(start,:), res(start,:)] = ...
        gmg(Mc, Lc, Ec, E_out_out, E_out_in, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, start, w, Eplot, uexact, MAX);
end

err_inf = err_inf(end:-1:1,:);
res = res(end:-1:1,:);

i = n_level-2;
sphplot = Eplot{i}*umg{i};
%sphplot = rhsfn(th_plot, phi_plot);
figure(5)
sphplot = reshape(sphplot, size(xp));
surf(xp, yp, zp, sphplot);
%title( '\fontsize{15} cos(\phi) ');
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
if n_level == 6
semilogy(n,res(1,:),'o--', n,res(2,:),'r*--', n,res(3,:),'g+--', n,res(4,:),'k-s', ...
         n,res(5,:),'c^-');
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
fs = 12;
set(gca,'Fontsize',fs)
%title('\fontsize{15} relative residuals in the \infty-norm')
xlabel('\fontsize{15} number of v-cycles')
ylabel('\fontsize{15} ||f^h-A^hu^h||_{\infty}/||f^h||_{\infty}')

% plot error of matlab and error of different number of vcycles
figure(3)
n = 1:MAX;
n = n-1;
if n_level == 6
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

% semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
% legend('N=20','N=10')
xlabel('\fontsize{15}number of vcyles')
% ylabel('\fontsize{15}|error|_{\infty}')
ylabel('\fontsize{15} ||u^h-u||_{\infty}/||u||_{\infty}')
%title('\fontsize{15}cos(3\theta)sin(\phi)^3(9cos(\phi)^2-1), Neumann B.C.s')
%title('\fontsize{15}cos(phi), Dirichlet B.C.s')
set(gca,'Fontsize',12)
