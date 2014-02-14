%% test geometric multigrid method to solve poisson equation on a quarter of a torus

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');


x0 = -3;
x1 = 3;

%%
% 2D example on a circle
% Construct a grid in the embedding space

%dx = 0.05; % grid size
dx = 0.0125;
dx_coarsest = 0.1;   % coarsest grid size
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

r = 0.5;
a = 1;
cen = [0,0,0];
cpf1 = @(x,y,z) cpQuarterTorus(x,y,z, a,r,cen);
cpf = @(x,y,z) cpbar_3d(x,y,z,cpf1);

has_boundary = true;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg] = ...
    build_mg_cpgrid3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);

n_level = length(a_band);

%disp('building cp matrices ... ')
%[Mc, Lc, Ec] = build_mg_cpmatrix3d(a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order);

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM3d(a_x1d, a_y1d, a_z1d, a_xcp, a_ycp, a_zcp, a_band, a_bdyg, p_f2c, p_c2f);


R = cell(n_level,1);
for i = 1:1:n_level
    R{i} = a_xg{i}.^2 + a_yg{i}.^2 + a_zg{i}.^2;
end


% pt is a point whose value we want to specify as the same at each level of
% V-Cycle.

pt = [1.4 0 0];

m = 1; n = 1;
% uexactfn = @(theta, phi) cos(4*m*theta).*exp(sin(n*phi));
% rhsfn = @(theta, phi) -16*m^2*cos(4*m*theta).*exp(sin(n*phi)) + ...
%          n^2*cos(4*m*theta).*exp(sin(n*phi)).*(cos(n*phi).^2-sin(n*phi));

uexactfn = @(theta, phi) cos(4*m*theta);
rhsfn = @(theta, phi, r) -16*m^2*cos(4*m*theta)./r.^2 - shift*uexactfn(theta,phi);

%% building E_plot for purpose of plotting and debug
% plotting grid on the torus 
[xp, yp, zp] = quarterTorus(r,50,a);

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


%% Debug setup
tol = 1e-10;
ind = cell(n_level,1);
for i = 1:1:n_level
    ind{i} = ( abs(a_xg{i}-pt(1))<tol & abs(a_yg{i}-pt(2))<tol & abs(a_zg{i}-pt(3))<tol );
end

[theta phi] = cart2paramTorus(xp,yp,zp,a);
uexact = uexactfn(theta(:),phi(:));
[th_pt phi_pt] = cart2paramTorus(pt(1),pt(2),pt(3),a);
uexact_pt = uexactfn(th_pt,phi_pt);

uexact_debug = cell(n_level,1);
for i = 1:1:n_level
    [theta phi] = cart2paramTorus(a_xcp{i}, a_ycp{i}, a_zcp{i},a);
    uexact_debug{i} = uexactfn(theta, phi);
end

%% Setting up right hand side
disp('building right hand side and allocate space for solution ... ')
F = cell(n_level,1);
V = cell(n_level,1);
for i = 1:1:n_level
    [theta phi] = cart2paramTorus(a_xcp{i}, a_ycp{i}, a_zcp{i},a);
    rad = a + r*cos(phi);
    F{i} = rhsfn(theta,phi,rad);
	V{i} = zeros(size(F{i}));
end
disp('done')

disp('making the problem to solve poisson equation on a circle well-posedness ... ')
tol = 1e-10;
Mc1 = Mc;
Lc1 = Lc;
F1 = F;
i = n_level;
%pt1 = pt + [dx_coarsest 2*dx_coarsest -dx_coarsest];
pt1 = pt;
j = abs(a_xg{i}-pt1(1)) < tol & abs(a_yg{i}-pt1(2)) < tol & abs(a_zg{i}-pt1(3)) < tol;
Lc1{i}(j,:) = sum(Ec{i});
Mc1{i}(j,:) = sum(Ec{i});
F1{i}(j) = 0;
for j = 1:1:length(a_xg{i})
            if abs(a_xg{i}(j)-pt1(1)) < tol && abs(a_yg{i}(j)-pt1(2)) < tol && abs(a_zg{i}(j)-pt1(3)) < tol
              Lc1{i}(j,:) = sum(Ec{i});
              Mc1{i}(j,:) = sum(Ec{i});
                disp(['changing the rows...  level:', num2str(i)]);
%                   Lc1{i}(j,:) = 0;
%                   %Lc1{i}(:,j) = 0;
%                   Lc1{i}(j,j) = 1;
%                   Mc1{i}(j,:) = 0;
%                   %Mc1{i}(:,j) = 0;
%                   Mc1{i}(j,j) = 1;
                  F1{i}(j) = 0;
            end
end

% Following line of code seems a bit confusing because a circle does not have a
% boundary, and the logical variable 'has_boundary' is set to be 'false',
% but in order to make the problem wellposedness, we do something with the
% coefficient matrix just as when dealing with Neumann Boundary Conditions,
% so we add the following line:
% [Mc, Lc, Ec, F] = app_bnd3d(Mc, Lc, Ec, F, a_xg, a_yg, a_zg, a_bdyg, pt, 'neumann');
% [Mc, Lc, Ec, F] = app_bnd3d(Mc, Lc, Ec, F, a_xg, a_yg, a_zg, a_bdyg, pt, 'dirichlet');
%disp('done')


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


u_matlab = cell(n_level,1);
if p == 3
    i = n_level;
    tic;
    unew = Mc{i} \ F{i};
    t_matlab1 = toc
    circplot = Eplot{i}*unew;
    error = circplot - uexact;
    error = error - (unew(ind{i}) - uexact_pt);
    err_inf_matlab1 = max(abs(error));
elseif p == 2
    i = n_level-1;
    tic;
    unew = Mc{i} \ F{i};
    t_matlab1 = toc
    circplot = Eplot{i}*unew;
    error = circplot - uexact;
    error = error - (unew(ind{i}) - uexact_pt);
    err_inf_matlab1 = max(abs(error));
end

%i = n_level-1;
%tic;
%unew = Mc{i} \ F{i};
%t_matlab2 = toc
%circplot = Eplot{i}*unew;
%err_inf_matlab2 = max(abs(uexact-circplot));


MAX = 50;
err_inf = zeros(n_level-1,MAX);   
res1 = zeros(n_level-1, MAX);
res2 = zeros(n_level-1, MAX);
err_matlab = zeros(n_level-1, MAX);
u_multigrid = cell(n_level-1,1);
u_mg_debug = cell(n_level-1,1);
for start = 1:1:n_level-1
    u_matlab{start} = zeros(size(a_xg{start}));
    tic; 
    [umg_cell err_inf(start,:) res1(start,:) res2(start,:) err_matlab(start,:) cnt] = ...
        gmg_debug(Mc1, Lc1, Ec, V, F1, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, u_matlab, Eplot, ind, uexact_pt, MAX);
%    [umg_cell err_inf(start,:) res1(start,:) res2(start,:) err_matlab(start,:) cnt] = ...
%        gmg_debug(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, u_matlab, Eplot, ind, uexact_pt, MAX);


    u_multigrid{start} = umg_cell{cnt};
    u_mg_debug{start} = umg_cell;

    t_mg = toc

end

%torusplot = Eplot{n_level-1}*umg;
torusplot = uexact;
figure(1)
torusplot = reshape(torusplot, size(xp));
surf(xp, yp, zp, torusplot);
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
if n_level == 4
semilogy(n,res1(1,:),'.-', n,res1(2,:),'r.-', n,res1(3,:),'c^-');
legend('N=40','N=20','N=10')
elseif n_level == 3
semilogy(n,res1(1,:),'.-', n,res1(2,:),'r.-');
legend('N=20','N=10')
end
% semilogy(n,res(1,:),'.-', n,res(2,:),'r.-');
% legend('N=20','N=10')
xlabel('number of vcyles')
ylabel('|residual|_{\infty}')
title(['sin(2\pi x) with p=', num2str(p), ',  res = E*(f-L*v)'])

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
if n_level == 4
semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-');
legend('N=40','N=20','N=10')
elseif n_level == 3
semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
legend('N=20','N=10')
end
% semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
% legend('N=20','N=10')
xlabel('number of vcyles')
ylabel('|error|_{\infty}')
title(['sin(2\pi x) with p=', num2str(p)])
