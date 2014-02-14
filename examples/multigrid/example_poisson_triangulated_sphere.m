% test geometric multigrid method to solve poisson equation on a triangulated 
% sphere, the right hand side function is only know at the vertices of the
% triangles, linear interpolation on each triangle is used to obtain the value
% of the right hand side function at the closest points lying somewhere in
% the triangles. The error is then measured on the true sphere.

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');

addpath('../../surfaces/tri')
addpath('../../surfaces/tri2cp')

% build the triangulated sphere using the code by Anton Semechko 
disp('build the triangulated sphere')
tic;
TR = IcosahedronMesh;
TR=SubdivideSphericalMesh(TR,8);
Faces = TR.Triangulation;
Vertices = TR.X;
toc
disp('done')


x0 = -3;
x1 = 3;

%%
% 2D example on a circle
% Construct a grid in the embedding space

%dx = 0.2; % grid size
dx = 0.025;
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

n_level = round(log(dx_coarsest/dx)/log(2)) + 1;
a_dx = cell(n_level,1);
a_x1d = cell(n_level,1);
a_y1d = cell(n_level,1);
a_z1d = cell(n_level,1);
a_band = cell(n_level,1);
a_xcp = cell(n_level,1);
a_ycp = cell(n_level,1);
a_zcp = cell(n_level,1);
a_dist = cell(n_level,1);
a_bdyg = cell(n_level,1);
a_xg = cell(n_level,1);
a_yg = cell(n_level,1);
a_zg = cell(n_level,1);
a_cp = cell(n_level,1);
a_cpface = cell(n_level,1);

disp('building cp grids ... ')
for i = 1:1:n_level
    a_dx{i} = dx*2^(i-1);
    a_x1d{i} = (x0:a_dx{i}:x1)';
    a_y1d{i} = a_x1d{i};
    a_z1d{i} = a_x1d{i};
    
    [ijk,dist,cp,xyz,cpface] = tri2cp(Faces, Vertices, a_dx{i}, x0);
    
    a_xcp{i} = cp(:,1);
    a_ycp{i} = cp(:,2);
    a_zcp{i} = cp(:,3);
    
    a_xg{i} = xyz(:,1);
    a_yg{i} = xyz(:,2);
    a_zg{i} = xyz(:,3);

    I = ijk(:,1);
    J = ijk(:,2);
    K = ijk(:,3);
    
    nx = length(a_x1d{i});
    ny = length(a_y1d{i});
    nz = length(a_z1d{i});
    a_band{i} = sub2ind([ny, nx, nz], J, I, K);
    
    a_dist{i} = dist;
    a_cp{i} = cp;
    a_cpface{i} = cpface;
end
disp('done')

disp('building cp matrices ... ')
Mc = cell(n_level,1);
Lc = cell(n_level,1);
Ec = cell(n_level,1);
for i = 1:1:n_level
    ddx = a_x1d{i}(2) - a_x1d{i}(1);
    Ec{i} = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, p);
    Ec{i} = Ec{i}(:, a_band{i});
    Lc{i} = laplacian_3d_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i});
    E = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, 1);
    E = E(:,a_band{i});
    Mc{i} = E*Lc{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});
end
disp('done')

shift = 1;
for i = 1:1:n_level
    Mc{i} = Mc{i} -shift*speye(size(Mc{i}));
    Lc{i} = Lc{i} -shift*speye(size(Lc{i}));
end

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM3d(a_x1d, a_y1d, a_z1d, a_xcp, a_ycp, a_zcp, a_band, a_bdyg, p_f2c, p_c2f);
disp('done')
 
% pt is a point whose value we want to specify as the same at each level of
% V-Cycle.

%pt = [1 0 0];
%uexactfn = @(th, phi) cos(phi+pi/2);
%rhsfn = @(th, phi) -2*cos(phi+pi/2) - shift*uexactfn(th, phi);

uexactfn = @(th, phi) cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1);
rhsfn = @(th, phi, r) ( -30*cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1) ) - shift*uexactfn(th,phi);

%% building E_plot for purpose of plotting and debug
% plotting grid on sphere, using theta as a parameterization
% cpxg = a_xcp{1};
% cpyg = a_ycp{1};
% cpzg = a_zcp{1};
% [th, phi, r] = cart2sph(cpxg,cpyg,cpzg);
% 
% u0 = uexactfn(th, phi);

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
F = cell(n_level,1);
V = cell(n_level,1);
[th, phi] = cart2sph(Vertices(:,1),Vertices(:,2),Vertices(:,3));
for i = 1:1:n_level
   f_vertices = rhsfn(th,phi);
   F{i} = interp1_tri(a_cp{i},a_cpface{i},Faces,Vertices,f_vertices);
   V{i} = zeros(size(F{i}));
end
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

i = n_level-2;
tic;
unew = Mc{i} \ F{i};
t_matlab3 = toc
circplot = Eplot{i}*unew;
err_inf_matlab3 = max(abs(uexact-circplot));

i = n_level-3;
tic;
unew = Mc{i} \ F{i};
t_matlab4 = toc
circplot = Eplot{i}*unew;
err_inf_matlab4 = max(abs(uexact-circplot));

% i = n_level-4;
% tic;
% unew = Mc{i} \ F{i};
% t_matlab5 = toc
% circplot = Eplot{i}*unew;
% err_inf_matlab5 = max(abs(uexact-circplot));

MAX = 50;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
R = [];
for start = 1:1:n_level-1
%     [thg rg] = cart2pol(a_xcp{start},a_ycp{start});
%     V{start} = v_initial_fn(thg);
%    F{start} = zeros(size(F{start}));
%    for i = start:1:n_level
%        V{i} = zeros(size(V{i}));
%    end
    tic;
    [umg err_inf(start,:) res(start,:)] = gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, Eplot, MAX);
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
if n_level == 5
semilogy(n,res(1,:),'.-', n,res(2,:),'r.-', n,res(3,:),'c^-', n,res(4,:),'k-s');
legend('N=40','N=20','N=10','N=5')
elseif n_level == 4
semilogy(n,res(1,:),'.-', n,res(2,:),'r.-', n,res(3,:),'c^-');
legend('N=20','N=10','N=5')
end
% semilogy(n,res(1,:),'.-', n,res(2,:),'r.-');
% legend('N=20','N=10')
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
if n_level == 5
semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-', n,err_inf(4,:), 'ko-');
legend('N=40','N=20','N=10','N=5')
elseif n_level == 4
semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-');
legend('N=20','N=10','N=5')
end
% semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
% legend('N=20','N=10')
xlabel('number of vcyles')
ylabel('|error|_{\infty}')
title(['cos(\phi) with p=', num2str(p)])

%how fine is the triangulation
for i = 100:1:200
vertices = Vertices(Faces(i,:),:);
edges = vertices(:,:)-vertices([2 3 1],:);
l_edges = sqrt(sum(edges.^2,2));
average = (sum(l_edges))/3
end

