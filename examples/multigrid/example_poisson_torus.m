%% test geometric multigrid method to solve poisson equation on a torus

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');

% add notay amg
addpath('/scratch/cheny1/opt/AGMG_3.1.1/Matlab')

x0 = -4;
x1 = 4;

%%
% 2D example on a circle
% Construct a grid in the embedding space

%dx = 0.05; % grid size
dx = 0.0125;
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

r = 0.6;
a = 1.2;
cen = [0,0,0];
cpf = @(x,y,z) cpTorus(x,y,z, a,r,cen);

has_boundary = false;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg] = ...
    build_mg_cpgrid3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);

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
   Mc{i} = E*Lc{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});
end
shift = 1;
for i = 1:1:n_level
    Mc{i} = Mc{i} - shift*speye(size(Mc{i}));
    Lc{i} = Lc{i} - shift*speye(size(Lc{i}));
end

%[Mc, Lc, Ec] = build_mg_cpmatrix3d(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order);

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM3d(a_x1d, a_y1d, a_z1d, a_xcp, a_ycp, a_zcp, a_band, a_bdyg, p_f2c, p_c2f);

 
R = cell(n_level,1);
for i = 1:1:n_level
    R{i} = a_xg{i}.^2 + a_yg{i}.^2 + a_zg{i}.^2;
end


% pt is a point whose value we want to specify as the same at each level of
% V-Cycle.

pt = [1.6 0 0];
%uexactfn = @(x,y,z) sin(2*pi*x);
%rhsfn = @(x,y,z) 

%uexactfn = @(x,y,z) z;
%rhsfn = @(x,y,z) 1/r^4./(x.^2+y.^2).*z.*( r^2*(-4*(x.^2+y.^2)+a*sqrt(x.^2+y.^2)) + 2*(x.^2+y.^2).*(a^2+x.^2+y.^2-2*a*sqrt(x.^2+y.^2)+z.^2) );

m = 3; k = 2;
uexactfn = @(theta, phi) sin(m*theta) + cos(k*phi);
rhsfn = @(theta, phi, rad) -m^2*sin(m*theta)./rad.^2 + ( -(-r*sin(phi))./(a+r*cos(phi))*k.*sin(k*phi) - k^2*cos(k*phi) )/r^2  - shift*uexactfn(theta, phi);

%k = 3;
%uexactfn = @(theta, phi) sin(k*phi);
%rhsfn = @(theta, phi, rad) ( (-r*sin(phi))./(a+r*cos(phi))*k.*cos(k*phi) - k^2*sin(k*phi) )/r^2  - shift*uexactfn(theta, phi);

%% building E_plot for purpose of plotting and debug
% plotting grid on the torus 
[xp, yp, zp] = torus(r,200,a);

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

%uexact = uexactfn(xp(:), yp(:), zp(:));
[theta,phi] = cart2paramTorus(xp(:), yp(:), zp(:), a);
uexact = uexactfn(theta,phi);

uexact_debug = cell(n_level,1);
for i = 1:1:n_level
    %uexact_debug{i} = uexactfn(a_xcp{i}, a_ycp{i}, a_zcp{i});
    [theta,phi] = cart2paramTorus(a_xcp{i}, a_ycp{i}, a_zcp{i},a);
    uexact_debug{i} = uexactfn(theta,phi);
end

disp('building right hand side and allocate space for solution ... ')
F = cell(n_level,1);
V = cell(n_level,1);
for i = 1:1:n_level
    %F{i} = rhsfn(a_xcp{i}, a_ycp{i}, a_zcp{i}); 	
	[theta, phi] = cart2paramTorus(a_xcp{i}, a_ycp{i}, a_zcp{i}, a);
    rad = a + r*cos(phi);
    F{i} = rhsfn(theta,phi,rad);
    V{i} = zeros(size(F{i}));
end
disp('done')

% disp('making the problem to solve poisson equation on a circle well-posedness ... ')
% tol = 1e-10;
% Mc1 = Mc;
% Lc1 = Lc;
% F1 = F;
% i = n_level;
% %pt1 = pt + [dx_coarsest 2*dx_coarsest -dx_coarsest];
% pt1 = pt;
% j = abs(a_xg{i}-pt1(1)) < tol & abs(a_yg{i}-pt1(2)) < tol & abs(a_zg{i}-pt1(3)) < tol;
% Lc1{i}(j,:) = sum(Ec{i});
% Mc1{i}(j,:) = sum(Ec{i});
% F1{i}(j) = 0;
% for j = 1:1:length(a_xg{i})
%             if abs(a_xg{i}(j)-pt1(1)) < tol && abs(a_yg{i}(j)-pt1(2)) < tol && abs(a_zg{i}(j)-pt1(3)) < tol
% %                 L{i}(j,:) = ones(1,length(L{i}));
% %                 M{i}(j,:) = ones(1,length(M{i}));
%                Lc1{i}(j,:) = sum(Ec{i});
%                Mc1{i}(j,:) = sum(Ec{i});
%                 disp(['changing the rows...  level:', num2str(i)]);
% %                  Lc1{i}(j,:) = 0;
% %                  Lc1{i}(:,j) = 0;
% %                  Lc1{i}(j,j) = 1;
% %                  Mc1{i}(j,:) = 0;
% %                  Mc1{i}(:,j) = 0;
% %                  Mc1{i}(j,j) = 1;
%                  F1{i}(j) = 0;
%             end
% end

% Following line of code seems a bit confusing because a circle does not have a
% boundary, and the logical variable 'has_boundary' is set to be 'false',
% but in order to make the problem wellposedness, we do something with the
% coefficient matrix just as when dealing with Neumann Boundary Conditions,
% so we add the following line:
% [Mc, Lc, Ec, F] = app_bnd3d(Mc, Lc, Ec, F, a_xcp, a_ycp, a_zcp, a_bdyg, 'neumann');
% [Mc, Lc, Ec, F] = app_bnd3d(Mc, Lc, Ec, F, a_xg, a_yg, a_zg, a_bdyg, pt, 'neumann');
% disp('done')


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

tol_agmg = 1e-6;

i = n_level;
tic;
%unew = Mc{i} \ F{i};
unew = agmg(Mc{i},F{i},[],tol_agmg,1000);
t_matlab0 = toc
circplot = Eplot{i}*unew;
err_inf_matlab0 = max(abs(uexact-circplot)) / norm(uexact,inf);

i = n_level-1;
tic;
%unew = Mc{i} \ F{i};
unew = agmg(Mc{i},F{i},[],tol_agmg,1000);
t_matlab1 = toc
circplot = Eplot{i}*unew;
err_inf_matlab1 = max(abs(uexact-circplot)) / norm(uexact,inf);

i = n_level-2;
tic;
%unew = Mc{i} \ F{i};
unew = agmg(Mc{i},F{i},[],tol_agmg,1000);
t_matlab2 = toc
circplot = Eplot{i}*unew;
err_inf_matlab2 = max(abs(uexact-circplot)) / norm(uexact,inf);

% i = n_level-3;
% tic;
% %unew = Mc{i} \ F{i};
% unew = agmg(Mc{i},F{i},[],tol_agmg,1000);
% t_matlab3 = toc
% circplot = Eplot{i}*unew;
% err_inf_matlab3 = max(abs(uexact-circplot)) / norm(uexact,inf);

MAX = 50;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
umg = cell(n_level,1);
for start = 1:1:n_level-1
    tic;
    [umg{start} err_inf(start,:) res(start,:)] = gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, Eplot, MAX);
    t_mg = toc
end
err_inf = err_inf(end:-1:1,:);
res = res(end:-1:1,:);

i = n_level-2;
torusplot = Eplot{i}*umg{i};
%torusplot = uexact;

h = figure(1);
torusplot = reshape(torusplot, size(xp));
handle = surf(xp, yp, zp, torusplot);
%title( '\fontsize{15} sin(3\theta)+cos(2\phi) ');
xlabel('\fontsize{15} x'); ylabel('\fontsize{15} y'); zlabel('\fontsize{15} z');
axis equal; shading interp;
colorbar;
fs = 15;
set(gca,'Fontsize',fs)
print(h,'-dpng','-r450','torus.png');

figure(2)
% rep_res_matlab = repmat(res_matlab, 1, 2);
% xx = [0 7];
% semilogy(xx,rep_res_matlab(1,:),'b',xx,rep_res_matlab(2,:),'r',xx,rep_res_matlab(3,:),'c', ...
%          xx,rep_res_matlab(4,:),'k');
% hold on

n = 0:MAX-1;
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

n = 0:MAX-1;
if n_level == 6
    semilogy(n,err_inf(1,:),'x-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-',n,err_inf(4,:),'k+-',n,err_inf(5,:),'m-d');
    legend('N=5','N=10','N=20','N=40','N=80')
elseif n_level == 5
    semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-',n,err_inf(4,:),'k+-');
    legend('N=10','N=20','N=40','N=80')
elseif n_level == 4
semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-');
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
%ylabel('\fontsize{15} |error|_{\infty}')
ylabel('\fontsize{15} ||u^h-u||_{\infty}/||u||_{\infty}')
%title('\fontsize{15} relative errors on a torus ')
