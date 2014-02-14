%% test geometric multigrid method to solve poisson equation on a semicircle

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');

x0 = -3;
x1 = 3;
y0 = -3;
y1 = 3;

%%
% 2D example on a semi-circle
% Construct a grid in the embedding space

dx = 0.003125/2; % grid size
dx_coarsest = 0.1;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';

dy = dx;

x1d = (x0:dx:x1)';
y1d = (y0:dx:y1)';

dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

n1 = 2;
n2 = 2;

p_f2c = 1;
p_c2f = 1;

w = 1;
radius = sqrt(3);
cpf1 = @(x,y) cpSemicircle(x,y,radius);  paramf = @paramSemicircle;  
%cpf = @(x,y) cpbar_2d(x,y,cpf1);
normal1 = [1,0];
normal2 = [1,0];
cpf = @(x,y) cptilde_openCurveIn2d(x,y,cpf1,normal1,normal2);

% If the curve or surface has 'real' boundary:
has_boundary = true;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);


disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg] = ...
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary, 0,0);

disp('building cp matrices ... ')
[Mc, Lc, Ec] = build_mg_cpmatrix(a_band, a_xcp, a_ycp, a_x1d, a_y1d, p, order);

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_band, a_bdyg, p_f2c, p_c2f);


n_level = length(a_band);
Eplot = cell(n_level-1,1);

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0, pi, 1000)';
r = radius*ones( size(thetas) );
% plotting grid in Cartesian coords
[xp, yp] = pol2cart(thetas, r);
xp = xp(:); yp = yp(:);
    
dx_tmp = dx;

for i = 1:1:n_level-1
    x = (x0:dx_tmp:x1)';
    y = (y0:dx_tmp:y1)';
    
    Eplot{i} = interp2_matrix_test( x, y, xp, yp, p );
    Eplot{i} = Eplot{i}(:, a_band{i});
    dx_tmp = 2*dx_tmp;
end
% for testing Dirichlet Boundary Conditions
shift = 0;
rhsfn = @(th,r) -2*ones(size(th))/radius^2;
uexactfn = @(th) th.*(pi-th); 

% k = 6;
% rhsfn = @(th,r)  k^2*exp( cos(k*th) ).*( sin(k*th).^2 - cos(k*th) )/radius^2;
% uexactfn = @(th) exp( cos(k*th) ) - exp(1);

% rhsfn = @(th,r) 4*exp( cos(2*th) ).*( sin(2*th).^2 - cos(2*th) );
% uexactfn = @(th) exp( cos(2*th) ) - exp(1);
% rhsfn = @(th,r) -sin(th);
% uexactfn = @(th) sin(th);
% rhsfn = @(th,r) -121*sin(11*th);
% uexactfn = @(th) sin(11*th);

% for testing Neumann Boundary Conditions
% shift = 1;
 %uexactfn = @(th) cos(10*th) - 1;
 %rhsfn = @(th,r) -100*cos(10*th) - shift*uexactfn(th);
 
% uexactfn = @(th) cos(th) - 1;
% rhsfn = @(th,r) -cos(th) - shift*uexactfn(th);
 
% rhsfn = @(th,r) -4*cos(2*th);
% uexactfn = @(th) cos(2*th) - 5;


uexact = uexactfn(thetas);

disp('building right hand side and allocate space for solution ... ')
[V, F] = helper_set_rhs(a_xcp, a_ycp, rhsfn, 1);
% for i = 1:1:n_level-1
%     V{i} = Ec{i}*V{i};
%     bdyg = logical(a_bdyg{i});
%     F{i}(bdyg,:) = -F{i}(bdyg,:);
% end
disp('done')

tol = 1e-10;
pt = [1 0];
%pt = [0.8 0.6];
%pt = [0.6 0.8];
%pt = [0 -1];

ind = cell(n_level,1);
for i = 1:1:n_level
    ind{i} = ( abs(a_xg{i}-pt(1))<tol & abs(a_yg{i}-pt(2))<tol );
end

[theta r] = cart2pol(pt(1),pt(2));
uexact_pt = uexactfn(theta);


disp('setting up boundary conditions ... ')
[Mc, Lc, Ec, F] = app_bnd(Mc, Lc, Ec, F, a_xcp, a_ycp, a_bdyg, pt, 'dirichlet');
%[Mc, Lc, Ec, F] = app_bnd(Mc, Lc, Ec, F, a_xg, a_yg, a_bdyg, pt,  'neumann');
for i = 1:1:n_level
    Mc{i} = Mc{i} - shift*speye(size(Mc{i}));
    Lc{i} = Lc{i} - shift*speye(size(Lc{i}));
end

disp('done')

semicircplot = cell(n_level-1,1);
error_inf_matlab = cell(n_level-1,1);
u_matlab = cell(n_level-1,1);
for i = 1:1:n_level-1
unew = Mc{i} \ F{i};
semicircplot{i} = Eplot{i}*unew;
error_inf_matlab{i} = max(abs( uexactfn(thetas) - semicircplot{i} )) / norm(uexactfn(thetas),inf);
u_matlab{i} = unew;
end


w = 1;
MAX = 50;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
res2 = zeros(n_level-1, MAX);
err_matlab = zeros(n_level-1, MAX);
u_multigrid = cell(n_level-1,1);
u_mg_debug = cell(n_level-1,1);

R = cell(n_level,1);

for start = 1:1:n_level-1
    [umg err_inf(start,:) res(start,:)] = ...
        gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, Eplot, MAX);
    %[umg_cell err_inf(start,:) res1(start,:) res2(start,:) err_matlab(start,:) cnt] = ...
    %    gmg_debug(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, u_matlab, Eplot, ind, uexact_pt, MAX);
    %u_multigrid{start} = umg_cell{cnt};
    %u_mg_debug{start} = umg_cell;
end

res1 = res;

% plot error of matlab and error of different number of vcycles
figure(1)

err_inf_matlab = cell2mat(error_inf_matlab);

rep_err_inf_matlab = repmat(err_inf_matlab,1,2);
xx = [0 7];
semilogy(xx,rep_err_inf_matlab(1,:),'b',xx,rep_err_inf_matlab(2,:),'r',xx,rep_err_inf_matlab(3,:),'c', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'g',xx,rep_err_inf_matlab(6,:),'m')
hold on

n = 0:MAX-1;
semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-',n,err_inf(4,:),'k-s',n,err_inf(5,:),'g+--',n,err_inf(6,:),'m-d');
legend('N=320','N=160','N=80','N=40','N=20','N=10')
%semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
%legend('N=20','N=10')
xlabel('number of vcyles')
ylabel('|error|_{\infty}')
% title(['cos(2*th)-5 with p=', num2str(p),  ';  Semicircle  Neumann  ; Relax by L, G-J ', num2str(n1), num2str(n2)])
% title(['cos(10*th)-101 with p=', num2str(p),  ';  Semicircle  Neumann  ; Relax by L, G-J ', num2str(n1), num2str(n2)])
% title(['sin(th) with p=', num2str(p),  ';  Semicircle  Dirichlet  ; Relax by L, G-J ', num2str(n1), num2str(n2)])
% title(['sin(11*th) with p=', num2str(p),  ';  Semicircle  Dirichlet  ; Relax by L, G-J ', num2str(n1), num2str(n2)])
% title(['exp(cos(2*th)) - exp(1) with p=', num2str(p), '; Semicircle  Dirichlet; Relax by L, G-J ', num2str(n1), num2str(n2)])

figure
if n_level == 9
    semilogy(n,res2(1,:),'.-',n,res2(2,:),'r.-',n,res2(3,:),'c^-', ...
             n,res2(4,:),'k-s',n,res2(5,:),'g+--',n,res2(6,:),'m-d', ...
             n,res2(7,:),'b*--', n,res2(8,:),'r*--');
    legend('N=1280', 'N = 640', 'N=320','N=160','N=80','N=40','N=20','N=10')
elseif n_level == 7
    semilogy(n,res2(1,:),'.-',n,res2(2,:),'r.-',n,res2(3,:),'c^-', ...
            n,res2(4,:),'k-s',n,res2(5,:),'g+--',n,res2(6,:),'m-d');
    legend('N=320','N=160','N=80','N=40','N=20','N=10')
    
end
% semilogy(n,res(1,:),'.-',n,res(2,:),'r*-');
% legend('N=20','N=10')
title('residual |f-M*u|')
xlabel('number of vcyles')
ylabel('|residual|_{\infty}')


figure
if n_level == 9
    semilogy(n,res1(1,:),'.-',n,res1(2,:),'r.-',n,res1(3,:),'c^-', ...
             n,res1(4,:),'k-s',n,res1(5,:),'g+--',n,res1(6,:),'m-d', ...
             n,res1(7,:),'b*--', n,res1(8,:),'r*--');
    legend('N=1280', 'N = 640', 'N=320','N=160','N=80','N=40','N=20','N=10')
elseif n_level == 7
    semilogy(n,res1(1,:),'.-',n,res1(2,:),'r.-',n,res1(3,:),'c^-', ...
            n,res1(4,:),'k-s',n,res1(5,:),'g+--',n,res1(6,:),'m-d');
    legend('N=320','N=160','N=80','N=40','N=20','N=10')
    
end
% semilogy(n,res(1,:),'.-',n,res(2,:),'r*-');
% legend('N=20','N=10')

title('residual |Eplot*(f-M*u)|')
xlabel('number of vcyles')
ylabel('|residual|_{\infty}')
%title(['sin(\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])
%title(['sin(\theta)+sin(',num2str(m),'\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])


% plot the grid of different level of the v-cycle
% figure(2);
% dx_tmp = dx;
% n_level = length(BAND);
% for i = 1:1:n_level
%    x = (-2.0:dx_tmp:2.0)';
%    y = x;
% 
%    nx = length(x);
%    ny = length(y);
% 
%    [xx yy] = meshgrid(x, y);
%    [cpx, cpy, dist, bdy] = cpbar_2d(xx, yy, cpf);
%    
%    bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
%    band = find(dist <= bw*dx_tmp);
%    
%    notband = setdiff(1:nx*ny, band);
%    
%    u = zeros(nx, ny);
%    [th, r] = cart2pol(xx, yy);
%    u = cos(th);
%    u = u(band);
%    
%    uplot = zeros(nx,ny);
%    uplot(band) = u;
%    uplot(notband) = NaN;
%    
%    subplot(3,2,i); hold off;
%    
%    % plot
%    pcolor(x,y,uplot);
% 
%    % make plot look pretty
%    caxis([-1.1 1.1]);
%    axis equal; colorbar;
%    title(['grid of the ' num2str(i) ' level of v-cycle'] );
%    xlabel('x'); ylabel('y');
%    hold on;
%    
%    tau = linspace(0,pi,100);
%    plot(cos(tau),sin(tau),'k');
%    
%    dx_tmp = dx_tmp*2;
% end
% 
% plotting grid on circle, using theta as a parameterization

