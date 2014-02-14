%% test geometric multigrid method to solve poisson equation on a Cos Curve
%% Parametrized by (x,cos(x)).

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');

% add notay amg
addpath('/scratch/cheny1/opt/AGMG_3.1.1/Matlab')

%%
% 2D example on a line segment
% Construct a grid in the embedding space

% relax by L.
% pt1 = [0.438744359656398   0.381558457093008];
% pt2 = [3.765516788149002   3.795199901137063];

 pt1 = [0.165648729499781   0.601981941401637];
 pt2 = [3.262971284540144   3.654079098476782];

% relax by M, after prolongation, do a closest point extension.

% not so good
% pt1 = [0.228976968716819   0.913337361501670];
% pt2 = [3.152378018969223   3.825816977489548];

% pt1 = [0.902716109915281   0.944787189721646];
% pt2 = [3.490864092468080   3.489252638400019];

% good
% pt1 = [0.817303220653433   0.868694705363510];
% pt2 = [3.084435845510910   3.399782649098896];

% pt1 = [0.853031117721894   0.622055131485066];
% pt2 = [3.350952380892271   3.513249539867053];

%pt1 = [rand rand];
%pt2 = [3+rand 3+rand];

% pt1 = [0,0];
% pt2 = [1,0];
% pt1(1) = pt1(1)-rand();
% pt2(1) = pt2(1)+rand();
a = pt1(1);
b = pt2(1);
kappa = (pt2(2)-pt1(2))/(pt2(1)-pt1(1));

normal1 = [pt2(2)-pt1(2) pt2(1)-pt1(1)];
normal1 = normal1/norm(normal1);
normal2 = normal1;

x0 = -2;
x1 = 6;
y0 = -2;
y1 = 6;

dx = 0.003125/2;
%dx = 0.00625;
dx_coarsest = 0.2;   % coarsest grid size
x1d_coarsest = x0 : dx_coarsest : x1;
y1d_coarsest = y0 : dx_coarsest : y1;


dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));


w = 1;

% Using "cpbar" [Macdonald, Brandman, Ruuth 2011]:
%cpf1 = @(x,y)cpCosCurve_d2c(x,y,endpt);
cpf1 = @(x,y)cpLineSegment(x,y,pt1,pt2);
cpf = @(x,y) cpbar_2d(x,y,cpf1);
%cpf = @(x,y) cptilde_openCurveIn2d(x,y,cpf1,normal1,normal2);

has_boundary = true;


disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg] = ...
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary,0,0);

%%
n1 = 3;
n2 = 3;

p_f2c = 1;
p_c2f = 1;

n_level = length(a_band);

disp('building cp matrices ... ')
[Mc, Lc, Ec] = build_mg_cpmatrix(a_band, a_xcp, a_ycp, a_x1d, a_y1d, p, order);

for i = 1:1:n_level
    ddx = a_x1d{i}(2) - a_x1d{i}(1);
    Ec{i} = interp2_matrix(a_x1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, p);
    Ec{i} = Ec{i}(:, a_band{i});
    Lc{i} = laplacian_2d_matrix(a_x1d{i}, a_y1d{i}, order, a_band{i}, a_band{i});
    E = interp2_matrix(a_x1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, 1);
    E = E(:,a_band{i});
    
    bdy = logical(a_bdyg{i});
    %E(bdy,:) = - E(bdy,:);
    Ec{i}(bdy,:) = - Ec{i}(bdy,:);
    %Lc{i}(bdy,:) = - Lc{i}(bdy,:);
     
%     [cpx_bar_double, cpy_bar_double] = cpbar_double_2d(a_xg{i}(bdy), a_yg{i}(bdy), cpf1);
%     E_bdy = interp2_matrix(a_x1d{i}, a_y1d{i}, cpx_bar_double, cpy_bar_double, p);
%     E_bdy = E_bdy(:, a_band{i});
%     Ec{i}(bdy,:) = E_bdy - 3*Ec{i}(bdy,:);
    
    Mc{i} = E*Lc{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});
    %I = speye(size(E));
    %Mc{i}(bdy,:) = (I(bdy,:)+Ec{i}(bdy,:));
    %Mc{i} = Ec{i}*Lc{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});

end

shift = 0;

for i = 1:1:n_level
    Mc{i} = Mc{i} - shift*speye(size(Mc{i}));
    Lc{i} = Lc{i} - shift*speye(size(Lc{i}));
end


disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_band, a_bdyg, p_f2c, p_c2f);



%% Neumann BC


%% Dirichlet BC
% Can not get the polynomial example work...
% k = 10;
% uexactfn = @(x) k*(x.^2 - (a+b)*x + a*b);
% rhsfn = @(x) 2*k*ones(size(x))/(1+kappa^2);

k1 = 7;
k = 2*k1*pi/(b-a);
uexactfn = @(x) exp( sin( k*(x-a) ) ) - 1;
rhsfn = @(x) k^2*exp(sin(k*(x-a))).*( cos(k*(x-a)).^2-sin(k*(x-a)) )/(1+kappa^2);

disp('building right hand side and allocate space for solution ... ')
F = cell(n_level,1);
V = cell(n_level,1);
for i = 1:1:n_level
    F{i} = rhsfn(a_xcp{i});
    bdy = logical(a_bdyg{i});
    %F{i}(bdy,:) = - F{i}(bdy,:);
    %F{i}(bdy,:) = 0;
    %V{i} = zeros(size(F{i}));
end
disp('done')


pt = [0, 1];
%[Mc, Lc, Ec, F] = app_bnd(Mc, Lc, Ec, F, a_xg, a_yg, a_bdyg, pt, 'neumann');
%[Mc, Lc, Ec, F] = app_bnd(Mc, Lc, Ec, F, a_xcp, a_ycp, a_bdyg, pt, 'dirichlet');

% for i = 1:1:n_level
%     bdy = logical(a_bdyg{i});
%     Lc{i}(bdy,:) = - Lc{i}(bdy,:);
% end

tol = 1e-10;
ind = cell(n_level,1);
for i = 1:1:n_level
    ind{i} = ( abs(a_xg{i}-pt(1))<tol & abs(a_yg{i}-pt(2))<tol );
end

uexact_pt = uexactfn(pt(1));

%% building E_plot for purpose of plotting and debug
%tt = linspace(a+dx_coarsest, b-dx_coarsest, 1000)';
tt = linspace(a, b, 1000)';
xp = tt;
yp = kappa*tt + (pt1(2)-kappa*pt1(1));

Eplot = cell(n_level-1,1);
dx_tmp = dx;
for i = 1:1:n_level-1
    x = (x0:dx_tmp:x1)';
    y = (y0:dx_tmp:y1)';
    
    Eplot{i} = interp2_matrix_test( x, y, xp, yp, p );
    Eplot{i} = Eplot{i}(:,a_band{i});
    dx_tmp = 2*dx_tmp;
end

uexact = uexactfn(tt);

circplot = cell(n_level-1,1);
error_inf_matlab = zeros(n_level-1,1);
res_matlab = zeros(n_level,1);
u_matlab = cell(n_level-1,1);
for i = 1:1:n_level-1
 unew = Mc{i} \ F{i};
 %unew = agmg(Mc{i},F{i},[],1e-12,1000);
 %unew = Ec{i}*unew;
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

end

matlab_order = log(error_inf_matlab(2:end)./error_inf_matlab(1:end-1))/log(2);
error_inf_matlab = error_inf_matlab(end:-1:1,:);
res_matlab = res_matlab(end:-1:1,:);
matlab_order = matlab_order(end:-1:1);

MAX = 100;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
u_multigrid = cell(n_level-1,1);
R = {};
for start = 1:1:n_level-1
    V{start} = zeros(size(F{start}));
    %V{start} = ones(size(F{start}));
    %V{start} = rand(size(F{start})) - 0.5;
    for i = start+1:1:n_level
        V{i} = zeros(size(F{i}));
    end
    [umg, err_inf(start,:), res(start,:)] = ...
        gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, start, w, uexact, Eplot, MAX);
    %[umg err_inf(start,:) res(start,:)] = ...
    %    gmg_debug(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, u_matlab, Eplot, ind, uexact_pt, MAX);
    u_multigrid{start} = umg;
end

err_inf = err_inf(end:-1:1,:);
res = res(end:-1:1,:);

figure(1)
% rep_res_matlab = repmat(res_matlab, 1, 2);
% xx = [0 7];
% semilogy(xx,rep_res_matlab(1,:),'b',xx,rep_res_matlab(2,:),'r',xx,rep_res_matlab(3,:),'c', ...
%          xx,rep_res_matlab(4,:),'k',xx,rep_res_matlab(5,:),'g',xx,rep_res_matlab(6,:),'m', ...
%          xx,rep_res_matlab(7,:),'--',xx,rep_res_matlab(8,:),'r--');
% hold on

n = 1:MAX;
n = n - 1;
if n_level == 8
    semilogy(n,res(1,:),'o--',n,res(2,:),'r*--',n,res(3,:),'g+--', ...
             n,res(4,:),'k-s',n,res(5,:),'c^-',n,res(6,:),'m-d', ...
             n,res(7,:),'b.-');
    legend('N=10','N=20','N=40','N=80','N=160','N=320','N=640')
elseif n_level == 6
    semilogy(n,res(1,:),'o--',n,res(2,:),'r*--',n,res(3,:),'g+--', ...
             n,res(4,:),'k-s',n,res(5,:),'c^-');
    legend('N=10','N=20','N=40','N=80','N=160')
elseif n_level == 4
    semilogy(n,res(1,:),'o--',n,res(2,:),'r*--',n,res(3,:),'g+--');
    legend('N=10','N=20','N=40')    
end
% semilogy(n,res(1,:),'.-',n,res(2,:),'r*-');
% legend('N=20','N=10')
fs = 12;
set(gca,'Fontsize',fs)
title('\fontsize{15} relative residuals in the \infty-norm')
xlabel('\fontsize{15} number of v-cycles')
ylabel('\fontsize{15} ||f^h-A^hu^h||_{\infty}/||f^h||_{\infty}')
%title('\fontsize{15} residual |Eplot*(f-A*u)|')
%xlabel('\fontsize{15} number of v-cycles')
%ylabel('\fontsize{15} |residual|_{\infty}')
%title(['sin(\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])
%title(['sin(\theta)+sin(',num2str(m),'\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])

% plot error of matlab and error of different number of vcycles
figure(2)

n = 1:MAX;
n = n - 1;
if n_level == 8
    semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', ...
         n,err_inf(4,:),'k-s',n,err_inf(5,:),'c^-',n,err_inf(6,:),'m-d', ...
            n,err_inf(7,:),'b.-');
    legend('N=10','N=20','N=40','N=80','N=160','N=320','N=640')
elseif n_level == 6
    semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', ...
         n,err_inf(4,:),'k-s',n,err_inf(5,:),'c^-');
    legend('N=10','N=20','N=40','N=80','N=160')
elseif n_level == 4
    semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--');
    legend('N=10','N=20','N=40')
end
hold on
%err_inf_matlab = cell2mat(error_inf_matlab);
rep_err_inf_matlab = repmat(error_inf_matlab,1,2);
xx = [0 7];
if n_level == 8
    semilogy(xx,rep_err_inf_matlab(1,:),'b--',xx,rep_err_inf_matlab(2,:),'r--',xx,rep_err_inf_matlab(3,:),'g', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'c', xx,rep_err_inf_matlab(6,:),'m-', ...
            xx,rep_err_inf_matlab(7,:),'b-');
elseif n_level == 6
    semilogy(xx,rep_err_inf_matlab(1,:),'b--',xx,rep_err_inf_matlab(2,:),'r--',xx,rep_err_inf_matlab(3,:),'g', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'c');
elseif n_level == 4
     semilogy(xx,rep_err_inf_matlab(1,:),'b--',xx,rep_err_inf_matlab(2,:),'r--',xx,rep_err_inf_matlab(3,:),'g');
end

% semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
% legend('N=20','N=10')

fs = 12;
set(gca,'Fontsize',fs)
xlabel('\fontsize{15} number of v-cycles')
ylabel('\fontsize{15} ||u^h-u||_{\infty}/||u||_{\infty}')
%title('\fontsize{14} relative errors on a cosine curve with Dirichlet B.C.s')
title('\fontsize{14} relative errors on a cosine curve with Neumann B.C.s')
%xlim([0 10])