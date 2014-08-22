%% test geometric multigrid method to solve poisson equation on a Bean like curve

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');


%%
% 2D example on a Bean shape Curve
% Construct a grid in the embedding space


x0 = -3;
x1 = 3;
y0 = -3;
y1 = 3;

%dx = 0.003125/2;
dx = 0.00625;
dx_coarsest = 0.2;   % coarsest grid size
x1d_coarsest = x0 : dx_coarsest : x1;
y1d_coarsest = y0 : dx_coarsest : y1;


dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));


w = 1; 

scale = 2;
pts = [...
         -0.4          0.42
         -0.6          0.4
         -0.8          0.25
         -0.9            0
         -0.8         -0.3
         -0.6         -0.43
         -0.3         -0.5
          0.1         -0.5
          0.4         -0.44
          0.7         -0.25
          0.8            0
          0.75          0.24
          0.6          0.35
          0.4          0.35
          0.17         0.25
         -0.05         0.25
         -0.22          0.35
         -0.4          0.42];
  
  % make it roughly centered at origin
  pts = pts + repmat([0.05 0.04], [size(pts,1), 1]);
  pts = scale*pts;

cpf = @(x,y) cpBeanCurve(x,y,scale);
has_boundary = false;
use_ndgrid = false;
need_param = true;

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg, a_param] = ...
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary, use_ndgrid, need_param);

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
    
    Mc{i} = E*Lc{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});
    %Mc{i} = Ec{i}*Lc{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});
    
    Mc_old{i} = Lc{i}*Ec{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});
    
%      [Ec{i}, GAMMA] = my_interp2_matrix(a_x1d{i},a_y1d{i},a_xcp{i},a_ycp{i},p,a_band{i});
%      Mc{i} = E*Lc{i} - GAMMA*(speye(size(E))-Ec{i});
end

% test whether M-matrix property holds
for i = 1:1:length(Mc)
diagM = diag(diag(Mc{i}));
flag1 = diagM > 0;
nnz(flag1)
tmp = Mc{i} - diagM;
flag2 = tmp < 0;
nnz(flag2)
max(abs(tmp(flag2)))
end

shift = 1;

for i = 1:1:n_level
    Mc{i} = Mc{i} - shift*speye(size(Mc{i}));
    Lc{i} = Lc{i} - shift*speye(size(Lc{i}));
end


disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_band, a_bdyg, p_f2c, p_c2f);

  
  sp = cscvn(pts');
  sp1 = fnder(sp);
  sp2 = fnder(sp,2);
  
  S1.type='()';
  S1.subs = {1,':'};
  S2.type='()';
  S2.subs = {2,':'};

  % parameterised curve:
  xs = @(t) subsref(ppval(sp,t), S1)';
  ys = @(t) subsref(ppval(sp,t), S2)';
  
  % derivative of parametrisation:
  xp = @(t) subsref(ppval(sp1,t), S1)';
  yp = @(t) subsref(ppval(sp1,t), S2)';
  
  % second derivative:
  xpp = @(t) subsref(ppval(sp2,t), S1)';
  ypp = @(t) subsref(ppval(sp2,t), S2)';


%% Setting up right hand side
[xtmp ytmp ttmp] = paramBeanCurve(200,scale);
total = max(ttmp);
coef = 2*pi/total;
%uexactfn = @(t) sin(coef*t);

%% ideal way is the following, which does not work on current matlab
%du_dt = matlabFunction( diff( sym(uexactfn) ) );
%du_ds = @(t) 1./sqrt(xp(t).^2+yp(t).^2).*du_dt(t);
%d2u_dsdt = matlabFunction( diff( sym(du_ds) ) );
%d2u_ds2 = @(t) 1./sqrt(xp(t).^2+yp(t).^2).*d2u_dsdt(t);

%% So we have to do the differentiation manually:
% %du_dt = @(t) coef*cos(coef*t);
% %du_ds = @(t) coef./sqrt(xp(t).^2+yp(t).^2).*cos(coef*t);

%d2u_dsdt = @(t) -coef^2*sin(coef*t)./sqrt(xp(t).^2+yp(t).^2) - coef*cos(coef*t).*(xp(t).*xpp(t)+yp(t).*ypp(t))./(xp(t).^2+yp(t).^2).^(1.5);
%d2u_ds2 = @(t) 1./sqrt(xp(t).^2+yp(t).^2).*d2u_dsdt(t);
%rhsfn = @(t) -shift*uexactfn(t) + d2u_ds2(t);

%% another exact solution and rhs function
uexactfn = @(t) sin(coef*t) + sin(10*coef*t);
rhsfn = @(t) -shift*uexactfn(t) - coef^2*(100*sin(10*coef*t)+sin(coef*t))./(xp(t).^2+yp(t).^2) - ...
                                  coef*(10*cos(10*coef*t)+cos(coef*t)).* (xp(t).*xpp(t)+yp(t).*ypp(t))./(xp(t).^2+yp(t).^2).^2;

disp('building right hand side and allocate space for solution ... ')
F = cell(n_level,1);
V = cell(n_level,1);
for i = 1:1:n_level
    F{i} = rhsfn(a_param{i});
    V{i} = zeros(size(F{i}));
end
disp('done')


pt = [0, 1];
%[Mc, Lc, Ec, F] = app_bnd(Mc, Lc, Ec, F, a_xg, a_yg, a_bdyg, pt, 'neumann');
%[Mc, Lc, Ec, F] = app_bnd(Mc, Lc, Ec, F, a_xcp, a_ycp, a_bdyg, pt, 'dirichlet');

tol = 1e-10;
ind = cell(n_level,1);
for i = 1:1:n_level
    ind{i} = ( abs(a_xg{i}-pt(1))<tol & abs(a_yg{i}-pt(2))<tol );
end

uexact_pt = uexactfn(pt(1));

%% building E_plot for purpose of plotting and debug
[xplot, yplot, tt] = paramBeanCurve(500,scale);

Eplot = cell(n_level-1,1);
dx_tmp = dx;
for i = 1:1:n_level-1
    x = (x0:dx_tmp:x1)';
    y = (y0:dx_tmp:y1)';
    
    Eplot{i} = interp2_matrix( x, y, xplot, yplot, p );
    Eplot{i} = Eplot{i}(:,a_band{i});
    dx_tmp = 2*dx_tmp;
end

uexact = uexactfn(tt');

circplot = cell(n_level-1,1);
error_inf_matlab = zeros(n_level-1,1);
res_matlab = zeros(n_level,1);
u_matlab = cell(n_level-1,1);
for i = 1:1:n_level-1
 unew = Mc{i} \ F{i};
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
error_inf_matlab = error_inf_matlab(end:-1:1);
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
        gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, Eplot, MAX);
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
            n,err_inf(7,:),'bx-');
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
xlim([0,10])
%title('\fontsize{15} relative errors on a bean-shaped curve')
%xlabel('\fontsize{15} number of v-cycles')
%ylabel('\fontsize{15} |error|_{\infty}')
%title('\fontsize{15} error on a bean-shaped curve ')
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
i = n_level-1;
figure(i);clf;plot2d_compdomain(u_multigrid{i},a_xg{i},a_yg{i},a_dx{i},a_dx{i},i)
axis([-2.3 2.3 -1.5 1.5])
[xp,yp] = paramBeanCurve(512,2);
hold on, plot(xp,yp,'k-','linewidth',2)
 set(gca,'Fontsize',12)
xlabel('\fontsize{15} x')
ylabel('\fontsize{15} y')
figure(5); plot(a_xg{i},a_yg{i},'r.','MarkerSize',10);
hold on, plot(a_xcp{i},a_ycp{i},'.');
set(gca,'Fontsize',12)
xlabel('\fontsize{15} x')
ylabel('\fontsize{15} y')
i = i-1; figure(6); plot(a_xg{i},a_yg{i},'r.','MarkerSize',7);
hold on, plot(a_xcp{i},a_ycp{i},'.');
set(gca,'Fontsize',12)
xlabel('\fontsize{15} x')
ylabel('\fontsize{15} y')



%% some tests:
for i = 1:n_level
Lc{i} = Lc{i} + shift*speye(size(Ec{i}));
Mc{i} = Mc{i} + shift*speye(size(Ec{i}));
end

for i = n_level:-1:1;
i
A = Ec{i} * (speye(size(Ec{i})) + 0.25*a_dx{i}^2*Lc{i});
for cnt = 1:10000
a = 2*(rand(length(A),1)-0.5);
%b = Ec{i}*a;
b = a;
%if norm(A*b,inf) > norm(b,inf)
if norm(A*b,inf) > norm(b,inf)    
disp('bad')
end
end
end