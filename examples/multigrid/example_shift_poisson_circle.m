%% test geometric multigrid method to solve -\Delta u + u = f on a circle

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');

% add notay amg
addpath('/scratch/cheny1/opt/AGMG_3.1.1/Matlab')

x0 = -3;
x1 = 3;
y0 = -3;
y1 = 3;

%%
% 2D example on a circle
% Construct a grid in the embedding space

dx = 0.003125;
%dx = 0.003125/2; % grid size
dx_coarsest = 0.2;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';

dy = dx;

dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
%bw = bw*2;

n1 = 3;
n2 = 3;

p_f2c = 1;
p_c2f = 1;

w = 1;

cpf = @cpCircle;

has_boundary = false;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg] = ...
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary, 0, 0);

n_level = length(a_band);

disp('building cp matrices ... ')
[Mc, Lc, Ec] = build_mg_cpmatrix(a_band, a_xcp, a_ycp, a_x1d, a_y1d, p, order);

w = 1-1e-5;
alpha = 0.69;  % seems that this is the critical weight that breaks the M-matrix property.
for i = 1:1:n_level
   ddx = a_x1d{i}(2) - a_x1d{i}(1);
   Ec{i} = interp2_matrix(a_x1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, p);
   %Ec{i} = interp2_matrix_cubicBspline(a_x1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, p);
   Ec{i} = Ec{i}(:, a_band{i});
   %Lc{i} = laplacian_2d_matrix(a_x1d{i}, a_y1d{i}, order, a_band{i}, a_band{i});
   Lc{i} = laplacian_wider_stencil_2d_matrix(a_x1d{i}, a_y1d{i}, order, alpha, 1-alpha, a_band{i}, a_band{i});
   E1 = interp2_matrix(a_x1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, 1,a_band{i});
   
   
   I = speye(size(Lc{i}));
   Mc{i} = E1*Lc{i} - 2*dim/ddx^2*(I-Ec{i});
   
%   my_E1 = my_interp2_matrix(a_x1d{i},a_y1d{i},a_xcp{i},a_ycp{i},nan, a_band{i});  
%   Mc{i} = E1*(Lc{i}+4/ddx^2*I) - 2*dim/ddx^2*I + 4/ddx^2*(Ec{i}-w*E1-(1-w)*my_E1);

   % What if we fix $\lambda$ at each grid level? Not converge.
   %if i < n_level
   %    Mc{i} = E*Lc{i} - 2*dim/a_dx{1}^2*(speye(size(Ec{i}))-Ec{i});
   %end
end

% test whether M-matrix property holds
for i = 1:1:length(Mc)
diagM = diag(diag(Mc{i}));
flag1 = diagM > 0;
nnz(flag1)
tmp = Mc{i} - diagM;
flag2 = tmp < 0;
nnz(flag2)

end

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_band, a_bdyg, p_f2c, p_c2f);


R = cell(n_level,1);
for i = 1:1:n_level
    R{i} = a_xg{i}.^2 + a_yg{i}.^2;
end
 
%% building E_plot for purpose of plotting and debug
% plotting grid on circle, using theta as a parameterization
thetas = linspace(0, 2*pi, 1000)';
r = ones( size(thetas) );
% plotting grid in Cartesian coords
[xp, yp] = pol2cart(thetas, r);
xp = xp(:); yp = yp(:);

Eplot = cell(n_level-1,1);
dxc = cell(n_level,1);
dx_tmp = dx;
for i = 1:1:n_level-1
    x = (x0:dx_tmp:x1)';
    y = (y0:dx_tmp:y1)';
    
    Eplot{i} = interp2_matrix( x, y, xp, yp, p );
    Eplot{i} = Eplot{i}(:,a_band{i});
    dxc{i} = dx_tmp;
    dx_tmp = 2*dx_tmp;
end
dxc{n_level} = dx_tmp;



%% Setting up right hand side
epsilon = 1;

%uexactfn = @(th,r) exp( cos(th) );
%rhsfn = @(th,r) exp( cos(th) ).*( sin(th).^2 - cos(th) ) - epsilon*uexactfn(th,r);

%uexactfn = @(th,r) exp( cos(11*th) ); 
%rhsfn = @(th,r) ( 121*exp( cos(11*th) ).*( sin(11*th).^2 - cos(11*th) ) ) - epsilon*uexactfn(th,r);

%uexactfn = @(th,r) exp( cos(2*th) );
%rhsfn = @(th,r) 4*exp( cos(2*th) ).*( sin(2*th).^2 - cos(2*th) ) - epsilon*uexactfn(th,r);

%uexactfn = @(th,r) log( 2 + cos(th) ); 
%rhsfn = @(th,r) ( -( 1 + 2*cos(th) ) ./ ( 2 + cos(th) ).^2 ) ./ (r.^2) - epsilon*uexactfn(th,r);

% n = 50;
% uexactfn = @(th,r) sin(th).^n;
% rhsfn = @(th,r) ( n*(n-1)*sin(th).^(n-2).*cos(th).^2 - n*sin(th).^n ) ./ (r.^2) - epsilon*uexactfn(th,r);

% n = 3;
% uexactfn = @(th,r) abs(sin(th)).^n;
% rhsfn = @(th,r) rhsfn_handle(th,r,n) - epsilon*uexactfn(th,r);

%n_mode = 50;
%uexactfn = @(th,r) uexactfn_mix_mode(th,n_mode);
%rhsfn = @(th,r) rhsfn_mix_mode(th,n_mode) - epsilon*uexactfn(th,r);

m = 12;
uexactfn = @(th,r) sin(th) + sin(m*th);
rhsfn = @(th,r) (-sin(th) - m^2*sin(m*th)) - epsilon*uexactfn(th,r);

  
  %uexactfn = @(th,r) sin(th) + sin(2*th) + sin(3*th) + sin(4*th) + sin(5*th);
  %rhsfn = @(th,r) -sin(th) - 4*sin(2*th) - 9*sin(3*th) - 16*sin(4*th) - 25*sin(5*th) - epsilon*uexactfn(th,r);

  %uexactfn = @(th,r) cos(th);
  %rhsfn = @(th,r) -cos(th) - epsilon*uexactfn(th,r);


%  uexactfn = @(th,r) sin(th);
%  rhsfn = @(th,r) -sin(th) - epsilon*uexactfn(th,r);

for i = 1:1:n_level
    Mc{i} = Mc{i} - epsilon*speye(size(Mc{i}));
    Lc{i} = Lc{i} - epsilon*speye(size(Lc{i}));
end

n_mode = 100;
random_coefficient = rand(n_mode,1);
uinitialfn = @(th) uinitial_mix_mode(th,n_mode,random_coefficient);

%uinitialfn = @(th) sin(th);

uexact = uexactfn(thetas);
%uexact = 0.5*ones(size(thetas));
%uexact = zeros(size(thetas));

uexact_debug = cell(n_level,1);
for i = 1:1:n_level
    [thg, rg] = cart2pol(a_xcp{i}, a_ycp{i});
    uexact_debug{i} = uexactfn(thg);
end
% v_initial_fn = @(th) sin(th);
% uexact = zeros(size(thetas));

tol = 1e-10;
%pt = [1 0];
pt = [0.8 0.6];
%pt = [0.6 0.8];
%pt = [0 -1];

ind = cell(n_level,1);
for i = 1:1:n_level
    ind{i} = ( abs(a_xg{i}-pt(1))<tol & abs(a_yg{i}-pt(2))<tol );
end

[theta r] = cart2pol(pt(1),pt(2));
uexact_pt = uexactfn(theta);


disp('building right hand side and allocate space for solution ... ')
[V, F] = helper_set_rhs(a_xcp, a_ycp, rhsfn, 1);
%[V, F] = helper_set_rhs(a_xg, a_yg, rhsfn, 1);
disp('done')


F_long = cell(n_level,1);
for i = 1:1:n_level
    F_long{i} = [F{i}; zeros(size(F{i}))];
end

Mc_rect = cell(n_level,1);
for i = 1:1:n_level
   lambda = 2*dim/dxc{i}^2;
   Mc_rect{i} = [Ec{i}*Lc{i}; lambda*(speye(size(Ec{i}))-Ec{i})];
end



N = cell(n_level,1);
N{1} = Mc{1};
N{n_level} = Mc{n_level};
for i = 1:1:n_level-2
    N{i+1} = TMf2c{i}*(Mc{i}*TMc2f{i});
end

circplot = cell(n_level-1,1);
error_inf_matlab = zeros(n_level-1,1);
res_matlab = zeros(n_level,1);
u_matlab = cell(n_level-1,1);
for i = 1:1:n_level-1
    tic;
    
    unew = Mc{i} \ F{i};
     
%    m = Mc{i};
%    p = amd(m);
%    x = m(p,p)\F{i}(p);
%    [Y,I] = sort(p);
%    unew = x(I);
    
    %unew = agmg(Mc{i},F{i},[],1e-12,1000);
    %unew = Mc_rect{i} \ F_long{i};
    
    %[L U] = ilu(Mc{i});
    %[unew flag relres iter] = bicgstab(Mc{i}, F{i}, 1e-10, 1000, L, U);
 
% By increasing the maximal number of iterations, bicgstab will converge;
% however as grid become finer, number of iterations will increase sharply
% [unew] = bicgstab(Mc{i}, F{i}, 1e-6, 1000);

% gmres seems not converge for this problem
% unew = gmres(Mc{i}, F{i}, 3, 1e-10, 200);

circplot{i} = Eplot{i}*unew;
error = circplot{i} - uexact;
%error = error - (unew(ind{i}) - uexact_pt );

error_inf_matlab(i) = max(abs( error )) / norm(uexact,inf);
res_matlab(i) = norm(Eplot{i}*(F{i} - Mc{i}*unew),inf);
u_matlab{i} = unew;

t_matlab = toc

end

matlab_order = log(error_inf_matlab(2:end)./error_inf_matlab(1:end-1))/log(2);
error_inf_matlab = error_inf_matlab(end:-1:1);
matlab_order = matlab_order(end:-1:1);

w = 1;
MAX = 40;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1,MAX);
res1 = zeros(n_level-1, MAX);
res2 = zeros(n_level-1, MAX);
err_matlab = zeros(n_level-1, MAX);
u_multigrid = cell(n_level-1,1);
u_mg_debug = cell(n_level-1,1);
for start = 1:1:n_level-1
    tic;
    [thg, rg] = cart2pol(a_xcp{start},a_ycp{start});
    %V{start} = uinitialfn(thg);
    %V{start} = uexactfn(thg);
    V{start} = zeros(size(F{start}));
    %V{start} = ones(size(F{start}));
    %V{start} = rand(size(F{start})) - 0.5;
    for i = start+1:1:n_level
        V{i} = zeros(size(F{i}));
    end
  
    [umg, err_inf(start,:), res1(start,:), res2(start,:), err_matlab(start,:), cnt] = ...
        gmg_debug(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, u_matlab, Eplot, ind, uexact_pt, MAX);
%    [umg, err_inf(start,:), res(start,:), n_vc, err_normal_v, err_normal_f] = ...
%        gmg_debug1(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, Eplot, MAX);
%    nn = length(err_normal_v(:));
    
    
%     figure 
%     subplot(2,1,1), semilogy(1:nn,err_normal_v(:),'.-','Linewidth',2,'Markersize',15);
%     hold on, semilogy(1:2*(n_level-start):nn, err_normal_v(1,:),'r.','Markersize',20);
%     xlabel('\fontsize{20} different round of V-Cycles and different grid levels')
%     ylabel('\fontsize{20} norm(v-E*v)')
%     fs = 15;
%     set(gca,'Fontsize',fs)
%     subplot(2,1,2), semilogy(1:nn,err_normal_f(:),'.-','Linewidth',2,'Markersize',15);
%     hold on, semilogy(1:2*(n_level-start):nn, err_normal_f(1,:),'r.','Markersize',20);
%     xlabel('\fontsize{20} different round of V-Cycles and different grid levels')
%     ylabel('\fontsize{20} norm(res-E*res)')
%     fs = 15;
%     set(gca,'Fontsize',fs)
%     
    
%     for i = 1:1:cnt
%         figure, plot(sqrt(a_xg{start}.^2+a_yg{start}.^2),umg{i}-uexact_debug{start},'.')
%         title(['finest dx: ', num2str(a_dx{start}), ';  ',num2str(i), '-th round of V-Cycle'])
%     end
    
    u_multigrid{start} = umg;
    u_mg_debug{start} = umg;
    
    t_mg = toc
end

err_inf = err_inf(end:-1:1,:);
res = res(end:-1:1,:);

res1 = res1(end:-1:1,:);

Mc_test = cell(n_level,1);
F_test = cell(n_level,1);
for i = 1:1:n_level
   Mc_test{i} = speye(size(Ec{i})) - Ec{i};
   flag = ( abs(a_xg{i}-a_xcp{i})+abs(a_yg{i}-a_ycp{i})<1e-10 );
   Mc_test{i}(flag,:) = Lc{i}(flag,:);
   F_test{i} = F{i};
   F_test{i}(logical(1-flag)) = 0;
end

figure
i = n_level-1;
dx = dxc{i};
tmp = Mc_test{i} \ F_test{i};
x1d = x0:dx:x1;
y1d = x1d;
[xx yy] = meshgrid(x1d,y1d);
zz = zeros(size(xx));
zz(a_band{i}) = tmp;
zz(a_band{i}) = uexact_debug{i};
%zz(a_band{i}) = F{i};
%zz(a_band{i}) = u_matlab{i};
%zz(a_band{i}) = u_multigrid{i};
contour(xx,yy,zz,32)
colorbar

% rep_res_matlab = repmat(res_matlab, 1, 2);
% xx = [0 7];
% semilogy(xx,rep_res_matlab(1,:),'b',xx,rep_res_matlab(2,:),'r',xx,rep_res_matlab(3,:),'c', ...
%          xx,rep_res_matlab(4,:),'k',xx,rep_res_matlab(5,:),'g',xx,rep_res_matlab(6,:),'m', ...
%          xx,rep_res_matlab(7,:),'--',xx,rep_res_matlab(8,:),'r--');
% hold on


n = 0:MAX-1;
figure
if n_level == 9
    semilogy(n,err_matlab(1,:),'o--',n,err_matlab(2,:),'r*--',n,err_matlab(3,:),'g+--', ...
             n,err_matlab(4,:),'k-s',n,err_matlab(5,:),'c^-',n,err_matlab(6,:),'m-d', ...
             n,err_matlab(7,:),'b.-', n,err_matlab(8,:),'r*-');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640','N=1280')
elseif n_level == 7
    semilogy(n,err_matlab(1,:),'o--',n,err_matlab(2,:),'r*--',n,err_matlab(3,:),'g+--', ...
            n,err_matlab(4,:),'k-s',n,err_matlab(5,:),'c^-',n,err_matlab(6,:),'m-d');
    legend('N=10','N=20','N=40','N=80','N=160','N=320')
elseif n_level == 8
    semilogy(n,err_matlab(1,:),'o--',n,err_matlab(2,:),'r*--',n,err_matlab(3,:),'g+--', ...
            n,err_matlab(4,:),'k-s',n,err_matlab(5,:),'c^-',n,err_matlab(6,:),'m-d', ...
            n,err_matlab(7,:),'b.-');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640')
end
title('algebraic error comparing to the solution of matlab backslash ')
xlabel('number of v-cycles')
ylabel('|u\_Vcycle-u\_Backslash|_{\infty}')

figure
if n_level == 9
    semilogy(n,res2(1,:),'--',n,res2(2,:),'r--',n,res2(3,:),'c^-', ...
             n,res2(4,:),'k-s',n,res2(5,:),'g-',n,res2(6,:),'m-d', ...
             n,res2(7,:),'b-', n,res2(8,:),'r-');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640','N=1280')
elseif n_level == 7
    semilogy(n,res2(1,:),'--',n,res2(2,:),'r--',n,res2(3,:),'c^-', ...
            n,res2(4,:),'k-s',n,res2(5,:),'g-',n,res2(6,:),'m-d');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320')
elseif n_level == 8
    semilogy(n,res2(1,:),'--',n,res2(2,:),'r--',n,res2(3,:),'c^-', ...
            n,res2(4,:),'k-s',n,res2(5,:),'g-',n,res2(6,:),'m-d', ...
            n,res2(7,:),'b-');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640')    
end
% semilogy(n,res(1,:),'.-',n,res(2,:),'r*-');
% legend('N=20','N=10')
title('relative residual in the \infty-norm')
xlabel('number of v-cycles')
ylabel('||f^h-A^hu^h||_{\infty}/||f^h||_{\infty}')


figure
if n_level == 9
    semilogy(n,res1(1,:),'o--',n,res1(2,:),'r*--',n,res1(3,:),'g+--', ...
             n,res1(4,:),'k-s',n,res1(5,:),'c^-',n,res1(6,:),'m-d', ...
             n,res1(7,:),'b.-', n,res1(8,:),'r.-');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640','N=1280')
elseif n_level == 7
    semilogy(n,res1(1,:),'o--',n,res1(2,:),'r*--',n,res1(3,:),'g+--', ...
            n,res1(4,:),'k-s',n,res1(5,:),'c^-',n,res1(6,:),'m-d');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320')
elseif n_level == 8
    semilogy(n,res1(1,:),'o--',n,res1(2,:),'r*--',n,res1(3,:),'g+--', ...
            n,res1(4,:),'k-s',n,res1(5,:),'c^-',n,res1(6,:),'m-d', ...
            n,res1(7,:),'bx-');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640')    
end
% semilogy(n,res(1,:),'.-',n,res(2,:),'r*-');
% legend('N=20','N=10')

fs = 12;
set(gca,'Fontsize',fs)

%title('\fontsize{15} relative residuals in the \infty-norm')
xlabel('\fontsize{15} number of v-cycles')
ylabel('\fontsize{15} ||f^h-A^hu^h||_{\infty}/||f^h||_{\infty}')
xlim([0,8])
%title('\fontsize{15} residual |Eplot*(f-A*u)|')
%xlabel('\fontsize{15} number of v-cycles')
%ylabel('\fontsize{15} |residual|_{\infty}')
%title(['sin(\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])
%title(['sin(\theta)+sin(',num2str(m),'\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])

% plot error of matlab and error of different number of vcycles
figure

%err_inf_matlab = cell2mat(error_inf_matlab);


if n_level == 9
    semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', ...
         n,err_inf(4,:),'k-s',n,err_inf(5,:),'c^-',n,err_inf(6,:),'m-d', ...
         n,err_inf(7,:),'b.-', n,err_inf(8,:),'r.-');
legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640','N=1280')
elseif n_level == 7
    semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', ...
        n,err_inf(4,:),'k-s',n,err_inf(5,:),'c^-',n,err_inf(6,:),'m-d');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320');
elseif n_level == 8
    semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', ...
        n,err_inf(4,:),'k-s',n,err_inf(5,:),'c^-',n,err_inf(6,:),'m-d',n,err_inf(7,:),'bx-');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640');
end
% semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
% legend('N=20','N=10')

hold on

rep_err_inf_matlab = repmat(error_inf_matlab,1,2);
xx = [0 6];
if n_level == 9
    semilogy(xx,rep_err_inf_matlab(1,:),'--',xx,rep_err_inf_matlab(2,:),'r--',xx,rep_err_inf_matlab(3,:),'g', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'c',xx,rep_err_inf_matlab(6,:),'m', ...
         xx,rep_err_inf_matlab(7,:),'b',xx,rep_err_inf_matlab(8,:),'r');
elseif n_level == 7
     semilogy(xx,rep_err_inf_matlab(1,:),'--',xx,rep_err_inf_matlab(2,:),'r--',xx,rep_err_inf_matlab(3,:),'g', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'c',xx,rep_err_inf_matlab(6,:),'m');
elseif n_level == 8
     semilogy(xx,rep_err_inf_matlab(1,:),'--',xx,rep_err_inf_matlab(2,:),'r',xx,rep_err_inf_matlab(3,:),'g', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'c',xx,rep_err_inf_matlab(6,:),'m', ...
         xx,rep_err_inf_matlab(7,:),'b');
end

fs = 12;
set(gca,'Fontsize',fs)
xlabel('\fontsize{15} number of v-cycles')
ylabel('\fontsize{15} ||u^h-u||_{\infty}/||u||_{\infty}')
%title('\fontsize{15} relative errors in the \infty-norm')
xlim([0,8])
%title('\fontsize{15} error comparing to the exact solution ')
%title(['sin(\theta) + 0.01*sin(100\theta)   with   p=', num2str(p)]);
%title(['sum up sin with modes from 1 to ', num2str(n_mode), '  with   p=', num2str(p)]);
%title(['sin(\theta)^{10}   with   p=', num2str(p)]);
%title(['log(2+cos(\theta))   with   p=', num2str(p)]);
%title(['exp(cos(\theta))   with   p=', num2str(p)])
%title(['exp(cos(2\theta))   with   p=', num2str(p)])
%title(['sin(\theta)   with   p=', num2str(p)])
%title(['sin(\theta)+sin(11\theta)   with   p=', num2str(p)])

% figure(2); set(gcf, 'Position', [410 700 800 800]);
%  
% circplot = Eplot*unew;
% circplot_mg = Eplot*umg;
% error_circ_inf_mg = max(abs( uexact - circplot_mg' ))
% error_circ_inf = max(abs( uexact - circplot' ))
% 
% plot(thetas, circplot, thetas, circplot_mg, 'b.', thetas, uexactfn(thetas),'r--');
% legend('solution by matlab', 'multigrid', 'exact solution')
% title('comparison between the exact solution and multigrid solution');



% % plot the grid of different level of the v-cycle
% figure(3);
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
%    [cpx, cpy, dist] = cpCircle(xx,yy);
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
%    subplot(3,3,i); hold off;
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
%    tau = linspace(0,2*pi,100);
%    plot(cos(tau),sin(tau),'k');
% 
%    
%    dx_tmp = dx_tmp*2;
% end
% 
L = cell(n_level-1,1);
% for i = 3:1:n_level-1
%     tic;
%     I = speye(size(Mc{i}));
%     E = Ec{i};
%     A = Mc{i};
%     L{i} = I;
% %     for k = 1:1:n1
% %         L{i} = E*(I-tril(A)\A)*L{i} + (I-E);
% %     end
%     L{i} = I - tril(A)\A;
%     toc;
% end

% i = n_level-1;
% P0 = TMf2c{i};
% Ik = TMc2f{i};
% %P = Mc{i+1} \ (P0*Mc{i}*Ec{i});
% P = Mc{i+1} \ (P0*Mc{i});
% C = speye(size(Mc{i})) - Ik*P;
% 
% % make P, Ik to be square matrix
% P_bar = [P;zeros(size(P,2)-size(P,1),size(P,2))];
% Ik_bar = [Ik zeros(size(Ik,1),size(Ik,1)-size(Ik,2))];
% 
% tic;
% I = speye(size(Mc{i}));
% E = Ec{i};
% A = Mc{i};
% L{i} = I;
% %     for k = 1:1:n1
% %         L{i} = E*(I-tril(A)\A)*L{i} + (I-E);
% %     end
% L{i} = I - tril(A)\A;
% toc;    
% 
% m = L{i}^n2*(I-Ik*P)*L{i}^n1;
% for k = 1:1:10
%     u1 = u_matlab{i};
%     u2 = u_mg_debug{i}{k};
%     norm(m*(u1-u2),inf)/norm((u1-u2),inf)
% end
% 
% n = P*L{i}^n1;
% for k = 1:1:10
%     u1 = u_matlab{i};
%     u2 = u_mg_debug{i}{k};
%     norm(n*(u1-u2),inf)/norm((u1-u2),inf)
% end
% 
% % P0 = TMf2c{i};
% % Ik = TMc2f{i};
% % P = Mc{i+1} \ (P0*Mc{i});
% % C = speye(size(Mc{i})) - Ik*P;
% 
% 
