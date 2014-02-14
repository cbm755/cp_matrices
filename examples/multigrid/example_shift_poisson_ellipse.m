%% test geometric multigrid method to solve -\Delta u + u = f on an ellipse

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');

x0 = -6;
x1 = 6;
y0 = -3;
y1 = 3;

%%
% 2D example on a circle
% Construct a grid in the embedding space


dx = 0.00625;
dx_coarsest = 0.2;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';

dy = dx;

x1d = (x0:dx:x1)';
y1d = (y0:dx:y1)';

dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

n1 = 3;
n2 = 3;

p_f2c = 1;
p_c2f = 1;

w = 1;

a = 4;
b = 1;
cen = [0 0];
cpf = @(x,y) cpEllipse(x,y,a,b,cen);

has_boundary = false;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg] = ...
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary,0,0);

n_level = length(a_band);

disp('building cp matrices ... ')
[Mc, Lc, Ec] = build_mg_cpmatrix(a_band, a_xcp, a_ycp, a_x1d, a_y1d, p, order);

for i = 1:1:n_level
    Mc{i} = Mc{i} - speye(size(Mc{i}));
    Lc{i} = Lc{i} - speye(size(Lc{i}));
end

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_band, a_bdyg, p_f2c, p_c2f);



%% Setting up right hand side

% it is wierd to pass into the argument 'r' in the rhsfn without using it,
% this results from some testing codes in the poisson equation on a circle, 
% TODO: fix this. 
uexactfn = @(t,r) sin(t);
rhsfn = @(t,r) -sin(t)./(a^2*sin(t).^2 + b^2*cos(t).^2) - (a^2-b^2)*sin(t).*cos(t).^2./(a^2*sin(t).^2 + b^2*cos(t).^2).^2 - uexactfn(t,r);

% n = 10;
% m = 1;
% uexactfn = @(t,r) sin(n*t) + sin(m*t);
% rhsfn = @(t,r) -n^2*sin(n*t)./(a^2*sin(t).^2 + b^2*cos(t).^2) - n*(a^2-b^2)*sin(t).*cos(t).*cos(n*t)./(a^2*sin(t).^2 + b^2*cos(t).^2).^2 ...
%                -m^2*sin(m*t)./(a^2*sin(t).^2 + b^2*cos(t).^2) - m*(a^2-b^2)*sin(t).*cos(t).*cos(m*t)./(a^2*sin(t).^2 + b^2*cos(t).^2).^2 ...
%                - uexactfn(t,r);

% v_initial_fn = @(th) sin(th);
% uexact = zeros(size(thetas));

tol = 1e-10;
pt = [a 0];
%pt = [0 b];
%pt = [0 -1];
%pt = [1 dx_coarsest];

ind = cell(n_level,1);
for i = 1:1:n_level
    ind{i} = ( abs(a_xg{i}-pt(1))<tol & abs(a_yg{i}-pt(2))<tol );
end

t = 0;
uexact_pt = uexactfn(t);

transformed_xcp = cell(n_level,1);
transformed_ycp = cell(n_level,1);
for i = 1:1:n_level
    transformed_xcp{i} = a_xcp{i}/a;
    transformed_ycp{i} = a_ycp{i}/b;
end

disp('building right hand side and allocate space for solution ... ')
[V, F] = helper_set_rhs(transformed_xcp, transformed_ycp, rhsfn, 1);
%[V, F] = helper_set_rhs(a_xg, a_yg, rhsfn, 1);
disp('done')


% the following code is only used for passing a parameter into the gmg
% function; it's due to some testing code in the circle case. hope to fix
% this..
R = cell(n_level,1);
for i = 1:1:n_level
    R{i} = a_xg{i}.^2 + a_yg{i}.^2;
end


%% building E_plot for purpose of plotting and debug
% plotting grid on circle, using theta as a parameterization
tt = linspace(0, 2*pi, 1000)';
xp = a*cos(tt);
yp = b*sin(tt);

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
error = error - (unew(ind{i}) - uexact_pt );

error_inf_matlab(i) = max(abs( error ));
res_matlab(i) = norm(Eplot{i}*(F{i} - Mc{i}*unew),inf);
u_matlab{i} = unew;

end

MAX = 50;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
u_multigrid = cell(n_level-1,1);
for start = 1:1:n_level-1
    [thg rg] = cart2pol(a_xcp{start},a_ycp{start});
    %V{start} = uinitialfn(thg);
    %V{start} = uexactfn(thg);
    V{start} = zeros(size(F{start}));
    %V{start} = ones(size(F{start}));
    %V{start} = rand(size(F{start})) - 0.5;
    for i = start+1:1:n_level
        V{i} = zeros(size(F{i}));
    end
    [umg err_inf(start,:) res(start,:)] = ...
        gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, Eplot, MAX);
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



% plot the grid of different level of the v-cycle
% figure;
% dx_tmp = dx_coarsest;
% n_level = length(a_band);
% k = 0;
% for i = n_level:-1:n_level-3
%    k = k+1;
%    
%    x = (x0:dx_tmp:x1)';
%    y = (y0:dx_tmp:y1)';
% 
%    nx = length(x);
%    ny = length(y);
% 
%    [xx yy] = meshgrid(x, y);
%    
%    notband = setdiff(1:nx*ny, a_band{i});
%    
%    u = zeros(ny, nx);
%    [th, r] = cart2pol(xx, yy);
%    u = cos(th);
%    u = u(a_band{i});
%    
%    uplot = zeros(ny,nx);
%    uplot(a_band{i}) = u;
%    uplot(notband) = NaN;
%    
%    subplot(2,2,k); hold off;
%    
%    % plot
%    pcolor(xx,yy,uplot);
% 
%    % make plot look pretty
%    caxis([-1.1 1.1]);
%    axis equal; colorbar;
%    title(['grid of the ' num2str(i) ' level of v-cycle'] );
%    xlabel('x'); ylabel('y');
%    hold on;
%    
%    tau = linspace(0,2*pi,100);
%    plot(a*cos(tau),b*sin(tau),'k');
% 
%    
%    dx_tmp = dx_tmp/2;
% end

