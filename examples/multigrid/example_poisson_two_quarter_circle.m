%% test geometric multigrid method to solve poisson equation on an 
%% open smooth curve which is two connected quarter of a circle; 
%% the shape is roughly like the Cos curve from 0 to pi.

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');


%%
% 2D example on a Cos Curve
% Construct a grid in the embedding space


x0 = -3;
x1 = 3;
y0 = -3;
y1 = 3;

dx = 0.003125;
dx_coarsest = 0.1;   % coarsest grid size
x1d_coarsest = (x0 : dx_coarsest : x1)';
y1d_coarsest = (y0 : dx_coarsest : y1)';


dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));


w = 1;

% Using "cpbar" [Macdonald, Brandman, Ruuth 2011]:
cpf1 = @(x,y)cpTwoQuarterCircles(x,y);
cpf = @(x,y) cpbar_2d(x,y,cpf1);

has_boundary = true;
 
disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg] = ...
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary,0,0);

n1 = 3;
n2 = 3;

p_f2c = 1;
p_c2f = 1;

n_level = length(a_band);

disp('building cp matrices ... ')
[Mc, Lc, Ec] = build_mg_cpmatrix(a_band, a_xcp, a_ycp, a_x1d, a_y1d, p, order);

% for i = 1:1:n_level
%     ddx = a_x1d{i}(2) - a_x1d{i}(1);
%     Ec{i} = interp2_matrix(a_x1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, p);
%     Ec{i} = Ec{i}(:, a_band{i});
%     Lc{i} = laplacian_2d_matrix(a_x1d{i}, a_y1d{i}, order, a_band{i}, a_band{i});
%     E = interp2_matrix(a_x1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, 1);
%     E = E(:,a_band{i});
%     
%     %bdy = logical(a_bdyg{i});
%     %E(bdy,:) = - E(bdy,:);
%     %Ec{i}(bdy,:) = - Ec{i}(bdy,:);
%     
%     %Mc{i} = E*Lc{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});
%     Mc{i} = Ec{i}*Lc{i} - 2*dim/ddx^2*(speye(size(E))-Ec{i});
% end

shift = 1;

for i = 1:1:n_level
    Mc{i} = Mc{i} - shift*speye(size(Mc{i}));
    Lc{i} = Lc{i} - shift*speye(size(Lc{i}));
end


disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_band, a_bdyg, p_f2c, p_c2f);



%% Setting up right hand side

%uexactfn = @(s) cos(s);
%rhsfn = @(s) -shift*uexactfn(s) - cos(s); 

%uexactfn = @(s) cos(s) + cos(2*s) + cos(3*s) + cos(4*s);
%rhsfn = @(s) -shift*uexactfn(s) - cos(s) - 4*cos(2*s) - 9*cos(3*s) - 16*cos(4*s); 

k = 10;
uexactfn = @(s) cos(s) + cos(k*s);
rhsfn = @(s) -shift*uexactfn(s) -cos(s) - k^2*cos(k*s); 

%uexactfn = @(s) sin(10*s);
%rhsfn = @(s) - 100*sin(10*s); 
%uexactfn = @(s) sin(s);
%rhsfn = @(s) - sin(s); 

disp('building right hand side and allocate space for solution ... ')
F = cell(n_level,1);
V = cell(n_level,1);
for i = 1:1:n_level
    s = zeros(size(a_xcp{i}));
    flag1 = (a_ycp{i}>=0);
    s(flag1) = pi/2-cart2pol(a_xcp{i}(flag1)+1,a_ycp{i}(flag1));
    flag2 = (a_ycp{i}<0);
    s(flag2) = cart2pol(a_xcp{i}(flag2)-1,a_ycp{i}(flag2))+3*pi/2;
    F{i} = rhsfn(s);
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
[xp yp tt] = paramTwoQuarterCircles(1000);

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
 %unew = Ec{i}*unew;
% By increasing the maximal number of iterations, bicgstab will converge;
% however as grid become finer, number of iterations will increase sharply
% [unew] = bicgstab(Mc{i}, F{i}, 1e-6, 1000);

% gmres seems not converge for this problem
% unew = gmres(Mc{i}, F{i}, 3, 1e-10, 200);

circplot{i} = Eplot{i}*unew;
error = circplot{i} - uexact;
%error = error - (unew(ind{i}) - uexact_pt );

error_inf_matlab(i) = max(abs( error ));
res_matlab(i) = norm(Eplot{i}*(F{i} - Mc{i}*unew),inf);
u_matlab{i} = unew;

end

MAX = 100;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
u_multigrid = cell(n_level-1,1);
R = {};
for start = 1:1:n_level-1
    V{start} = zeros(size(F{start}));
    %V{start} = u_matlab{start};
    %V{start} = ones(size(F{start}));
    %V{start} = rand(size(F{start})) - 0.5;
    for i = start+1:1:n_level
        V{i} = zeros(size(F{i}));
    end
    [umg err_inf(start,:) res(start,:)] = ...
        gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, Eplot, MAX);
    %[umg err_inf(start,:) res(start,:)] = ...
    %    gmg_debug(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, u_matlab, Eplot, ind, uexact_pt, MAX);
    u_multigrid{start} = umg;
end

figure(1)
% rep_res_matlab = repmat(res_matlab, 1, 2);
% xx = [0 7];
% semilogy(xx,rep_res_matlab(1,:),'b',xx,rep_res_matlab(2,:),'r',xx,rep_res_matlab(3,:),'c', ...
%          xx,rep_res_matlab(4,:),'k',xx,rep_res_matlab(5,:),'g',xx,rep_res_matlab(6,:),'m', ...
%          xx,rep_res_matlab(7,:),'--',xx,rep_res_matlab(8,:),'r--');
% hold on

n = 1:MAX;
n = n - 1;
if n_level == 6
    semilogy(n,res(1,:),'.-',n,res(2,:),'r.-',n,res(3,:),'c^-', ...
             n,res(4,:),'k-s',n,res(5,:),'g+--');
    legend('N=320','N=160','N=80','N=40','N=20')
elseif n_level == 4
    semilogy(n,res(1,:),'.-',n,res(2,:),'r.-',n,res(3,:),'c^-');
    legend('N=80','N=40','N=20')
    
end
% semilogy(n,res(1,:),'.-',n,res(2,:),'r*-');
% legend('N=20','N=10')

xlabel('number of vcyles')
ylabel('|residual|_{\infty}')
%title(['sin(\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])
%title(['sin(\theta)+sin(',num2str(m),'\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])

% plot error of matlab and error of different number of vcycles
figure(2)

%err_inf_matlab = cell2mat(error_inf_matlab);
rep_err_inf_matlab = repmat(error_inf_matlab,1,2);
xx = [0 7];
if n_level == 6
    semilogy(xx,rep_err_inf_matlab(1,:),'b',xx,rep_err_inf_matlab(2,:),'r',xx,rep_err_inf_matlab(3,:),'c', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'g');
elseif n_level == 4
     semilogy(xx,rep_err_inf_matlab(1,:),'b',xx,rep_err_inf_matlab(2,:),'r',xx,rep_err_inf_matlab(3,:),'c');
end
hold on

n = 1:MAX;
n = n - 1;
if n_level == 6
    semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-', ...
         n,err_inf(4,:),'k-s',n,err_inf(5,:),'g+--');
legend('N=320','N=160','N=80','N=40','N=20')
elseif n_level == 4
    semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-');
    legend('N=80','N=40','N=20');
end
% semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
% legend('N=20','N=10')

xlabel('number of vcyles')
ylabel('|error|_{\infty}')
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
