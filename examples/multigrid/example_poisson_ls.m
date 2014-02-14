%% test geometric multigrid method to solve poisson equation on a surface
%% defined by the level set function (x-z^2)^2+y^2+z^2-1

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');


phi = @(x,y,z) (x-z.^2).^2+y.^2+z.^2-1;

disp('Construct anonymous function needed in the CPop function written by Tom')

syms sx sy sz sh; 
sphi = phi( sx, sy, sz );
grad_sphi = jacobian( sphi,[sx sy sz] );

sf = -grad_sphi; 
sf = sf/(sf * transpose(sf));
sf = simplify(transpose(sf)) +  sh*[sx ; sy ; sz];  % ... +sh*[sx ; sy ; sz] tri

% back to matlab anonymous functions
f = matlabFunction(sf,'vars',[sx sy sz sh]);
disp('done')

disp('Constructing the rhs function from the exact solution')
uexactfn = @(x,y,z) x.*y;
rhsfn = laplace_beltrami_ls3d(phi,uexactfn);
disp('done')

x0 = -2.2;
x1 = 2.2;

%%
% 2D example on a circle
% Construct a grid in the embedding space

%dx = 0.1; % grid size
dx = 0.0125;   % roughly OK.
%dx = 0.00625; % takes half an hour to compute CPs for the finest level
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

cpf = @(x,y,z) cpLevelSet3d(x,y,z,phi,f);
cpf1 = @(x,y,z) cp_ls(x,y,z,cpf);
cpf2 = @(x,y,z) cp_ls_test(x,y,z,cpf,50);

has_boundary = false;

%matlabpool(4)
matlabpool(12)

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg] = ...
    build_mg_cpgrid3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf2, has_boundary);
disp('done');

matlabpool close

n_level = length(a_band);

disp('building cp matrices ... ')
Mc = cell(n_level,1);
Lc = cell(n_level,1);
Ec = cell(n_level,1);
for i = 1:1:n_level
   ddx = a_x1d{i}(2) - a_x1d{i}(1);
   Ec{i} = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, p);
   Ec{i} = Ec{i}(:, a_band{i});
   %[Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_band{i});
   %E = interp3_matrix(a_x1d{i}, a_y1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, 2);
   %E = E(:,a_band{i});
   %Lc{i} = Dxc*(E*Dxc) + Dyc*(E*Dyc) + Dzc*(E*Dzc);
   Lc{i} = laplacian_3d_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i});
   E = interp3_matrix(a_x1d{i}, a_y1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, 1);
   E = E(:,a_band{i});
   Mc{i} = E*Lc{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});
end
shift = 1;
for i = 1:1:n_level
    Mc{i} = Mc{i} - shift*speye(size(Mc{i}));
    Lc{i} = Lc{i} - shift*speye(size(Lc{i}));
end
disp('done');
rhsfn = @(x,y,z) rhsfn(x,y,z) - shift*uexactfn(x,y,z);

%[Mc, Lc, Ec] = build_mg_cpmatrix3d(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order);

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM3d(a_x1d, a_y1d, a_z1d, a_xcp, a_ycp, a_zcp, a_band, a_bdyg, p_f2c, p_c2f);
disp('done')


%% building E_plot for purpose of plotting and debug
% plotting grid for the level set surface 
[xp yp zp] = paramLevelSet1(200);

n_level = length(a_band);
Eplot = cell(n_level,1);
dx_tmp = dx;
for i = 1:1:n_level
    x = (x0:dx_tmp:x1)';
    y = x;
    z = x;
    
    Eplot{i} = interp3_matrix( x, y, z, xp(:), yp(:), zp(:), p );
    Eplot{i} = Eplot{i}(:,a_band{i});
    dx_tmp = 2*dx_tmp;
end

%% Setting up the exact solution.

uexact = uexactfn(xp,yp,zp);

uexact_debug = cell(n_level,1);
for i = 1:1:n_level
    uexact_debug{i} = uexactfn(a_xcp{i},a_ycp{i},a_zcp{i});
end

disp('building right hand side and allocate space for solution ... ')
F = cell(n_level,1);
V = cell(n_level,1);
for i = 1:1:n_level
    F{i} = rhsfn(a_xcp{i},a_ycp{i},a_zcp{i});
    V{i} = zeros(size(F{i}));
end
disp('done')



i = n_level;
tic;
unew = Mc{i} \ F{i};
t_matlab0 = toc
circplot = Eplot{i}*unew;
err_inf_matlab0 = max(abs(uexact(:)-circplot));
err_inf_matlab0 = err_inf_matlab0/max(abs(uexact(:)));

i = n_level-1;
tic;
unew = Mc{i} \ F{i};
t_matlab1 = toc
circplot = Eplot{i}*unew;
err_inf_matlab1 = max(abs(uexact(:)-circplot));
err_inf_matlab1 = err_inf_matlab1/max(abs(uexact(:)));

i = n_level-2;
tic;
unew = Mc{i} \ F{i};
t_matlab2 = toc
circplot = Eplot{i}*unew;
err_inf_matlab2 = max(abs(uexact(:)-circplot));
err_inf_matlab2 = err_inf_matlab2/max(abs(uexact(:)));

%%
% 1 min
% err = 6.572344764553499e-04
% i = n_level-3;
% tic;
% unew = Mc{i} \ F{i};
% t_matlab3 = toc
% circplot = Eplot{i}*unew;
% err_inf_matlab3 = max(abs(uexact(:)-circplot));
% err_inf_matlab3 = err_inf_matlab3/max(abs(uexact(:)));

MAX = 10;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
DMc = cell(n_level,1);
tic;
for i = 1:1:n_level-1
    DMc{i} = diag(Mc{i});
end
toc;
umg = cell(n_level,1);
for start = 1:1:n_level-1
    tic;
    %[umg err_inf(start,:) res(start,:)] = gmg_t(Mc, DMc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, [], n1, n2, start, w, uexact(:), Eplot, MAX);
    [umg{start} err_inf(start,:) res(start,:)] = gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, [], n1, n2, start, w, uexact(:), Eplot, MAX);
    t_mg = toc
end
res = res(end:-1:1,:);
err_inf = err_inf(end:-1:1,:);


i = n_level-2;
lsplot = Eplot{i}*umg{i};
%torusplot = uexact;
figure1 = figure(1);
lsplot = reshape(lsplot, size(xp));
handle = surf(xp, yp, zp, lsplot);
%title( '\fontsize{18} u(x,y,z)=xy on \phi=(x-z^2)^2+y^2+z^2-1=0');
xlabel('\fontsize{18} x'); ylabel('\fontsize{18} y'); zlabel('\fontsize{18} z');
axis equal; shading interp; %camlight left
colorbar;
fs = 18;
set(gca,'Fontsize',fs)
%print(figure1,'-dpng','-r450','ls-dziuk-camlight')
print(figure1,'-dpng','-r450','ls-dziuk')

figure(2)
% rep_res_matlab = repmat(res_matlab, 1, 2);
% xx = [0 7];
% semilogy(xx,rep_res_matlab(1,:),'b',xx,rep_res_matlab(2,:),'r',xx,rep_res_matlab(3,:),'c', ...
%          xx,rep_res_matlab(4,:),'k');
% hold on

n = 1:MAX;
n = n-1;
if n_level == 5
semilogy(n,res(1,:),'.-', n,res(2,:),'r.-', n,res(3,:),'c^-', n,res(4,:),'k+-');
legend('N=10','N=20','N=40','N=80')
elseif n_level == 4
semilogy(n,res(1,:),'.-', n,res(2,:),'r.-', n,res(3,:),'c^-');
legend('N=10','N=20','N=40')
elseif n_level == 3
semilogy(n,res(1,:),'.-', n,res(2,:),'r.-');
legend('N=10','N=20')
end
% semilogy(n,res(1,:),'.-', n,res(2,:),'r.-');
% legend('N=20','N=10')
fs = 12;
set(gca,'Fontsize',fs)

title('\fontsize{15} residual |Eplot*(f-A*u)|')
xlabel('\fontsize{15} number of v-cycles')
ylabel('\fontsize{15} |residual|_{\infty}')

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
n = n-1;
if n_level == 6
    semilogy(n,err_inf(1,:),'x-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-',n,err_inf(4,:),'k+-',n,err_inf(5,:),'m-d');
    legend('N=10','N=20','N=40','N=80','N=160')
elseif n_level == 5
    semilogy(n,err_inf(1,:),'x-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-',n,err_inf(4,:),'k+-');
    legend('N=10','N=20','N=40','N=80')
elseif n_level == 4
semilogy(n,err_inf(1,:),'x-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-');
legend('N=10','N=20','N=40')
elseif n_level == 3
semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
legend('N=10','N=20')
end
% semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
% legend('N=20','N=10')
fs = 12;
set(gca,'Fontsize',fs)
xlabel('\fontsize{15} number of v-cycles')
ylabel('\fontsize{15} ||u^h-u||_{\infty}/||u||_{\infty}')
xlim([0,10])
ylim([5e-5,1])
%title('\fontsize{15} relative errors on \phi=(x-z^2)^2+y^2+z^2-1=0')
