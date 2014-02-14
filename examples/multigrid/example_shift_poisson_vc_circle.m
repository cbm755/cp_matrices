%% test geometric multigrid method to solve poisson equation on a circle

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

dx = 0.003125; % grid size
dx_coarsest = 0.2;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';

dy = dx;

x1d = (x0:dx:x1)';
y1d = (y0:dx:y1)';

dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: here we only implement the second order scheme
            % for the variable coefficient case.

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

n1 = 3;
n2 = 3;

p_f2c = 1;
p_c2f = 1;

w = 1;

cpf = @cpCircle;

has_boundary = false;

% beta_fun = @(theta,r) ones(size(theta));
% rhsfn = @(th,r) exp( cos(th) ).*( sin(th).^2 - cos(th) );
% uexactfn = @(th) exp( cos(th) );
% rhsfn = @(th,r) ( 121*exp( cos(11*th) ).*( sin(11*th).^2 - cos(11*th) ) );
% uexactfn = @(th) exp( cos(11*th) );  

 C = 1.5;
 n = 8; m = 5;
 beta_fun = @(theta,r) sin(n*theta) + C;
 rhsfn = @(th,r) n*cos(n*th).*cos(th) - sin(n*th).*sin(th) + n*cos(m*th).*cos(n*th) - m*sin(n*th).*sin(m*th) - C*sin(th) - C*m*sin(m*th);
 uexactfn = @(th) sin(th) + sin(m*th)/m;

%n = 5; m = 3;
%beta_fun = @(theta,r) exp(sin(n*theta));
%uexactfn = @(th) exp(sin(m*th));
%rhsfn = @(theta,r) m*exp(sin(m*theta)+sin(n*theta)).*(m*cos(m*theta).^2-m*sin(m*theta)+n*cos(m*theta).*cos(n*theta));

epsilon = 1;
rhsfn = @(theta,r) rhsfn(theta,r) - epsilon*uexactfn(theta);

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg] = ...
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary, 0, 0);

n_level = length(a_band);

disp('building cp matrices ... ')
Ec = cell(n_level,1);
Lc = cell(n_level,1);
Mc = cell(n_level,1);

for i = 1:1:n_level
	
	x = a_x1d{i};
	y = a_y1d{i};
	dx1 = x(2) - x(1);
	dy1 = y(2) - y(1);
    E1 = interp2_matrix(x,y,a_xcp{i},a_ycp{i},1);
    E1 = E1(:,a_band{i});
    E2 = interp2_matrix(x,y,a_xcp{i},a_ycp{i},2);
    E2 = E2(:,a_band{i});
	Ec{i} = interp2_matrix(x, y, a_xcp{i}, a_ycp{i}, p);
	Ec{i} = Ec{i}(:,a_band{i});

    % If the coefficients were constant, just call 'laplacian_2d_matrix' then done;
	% for the variable coefficient case, we need a few more work to do...
    
    % First compute the forward and backward difference in the x and y
    % direction
    [Dxb, Dxf, Dyb, Dyf] = firstderiv_upw1_2d_matrices(x,y,a_band{i});
    
    % Then compute the coefficients: beta_{i+1/2,j}, beta_{i-1/2,j},
    % beta_{i,j+1/2}, and beta_{i,j-1/2}.
    [theta, r] = cart2pol(a_xcp{i}, a_ycp{i});
    beta = beta_fun(theta, r);
    %[Axb, Axf, Ayb, Ayf] = sum_nb_2d_matrices(x,y,a_band{i});
    [Axb, Axf, Ayb, Ayf] = average_upw1_2d_matrices(x,y,a_band{i});
    coef_xb = Axb*beta; 
    coef_xf = Axf*beta;
    coef_yb = Ayb*beta;
    coef_yf = Ayf*beta;
    % We just want to multiply the value of beta with each corresponding row of 
    % the difference matrices, but making the vector beta into a sparse
    % diagonal matrix, and then multiply the difference matrices might be
    % faster; another choice would be 'bsxfun'; better way to do this?
	s = length(a_band{i});
    diag_coef_xb = spdiags(coef_xb,0,s,s); 
    diag_coef_xf = spdiags(coef_xf,0,s,s);
    diag_coef_yb = spdiags(coef_yb,0,s,s);
    diag_coef_yf = spdiags(coef_yf,0,s,s);
	
	Lc{i} = (diag_coef_xf*Dxf - diag_coef_xb*Dxb)/dx1 + ...
	        (diag_coef_yf*Dyf - diag_coef_yb*Dyb)/dy1;

    I = speye(size(Lc{i}));
    
    % E3*L is better than E1*L
    Mc{i} = E1*Lc{i}- 2*dim/a_dx{i}^2*(I-Ec{i});
	
    % Because of the shift, E3*L*E3 also works
    %Mc{i} = Ec{i}*Lc{i}*Ec{i};
    
    Lc{i} = Lc{i} - epsilon*I;
    Mc{i} = Mc{i} - epsilon*I;
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
    
    Eplot{i} = interp2_matrix_test( x, y, xp, yp, p );
    Eplot{i} = Eplot{i}(:,a_band{i});
    dxc{i} = dx_tmp;
    dx_tmp = 2*dx_tmp;
end
dxc{n_level} = dx_tmp;

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
    
    %m = Mc{i};
    %p = amd(m);
    %x = m(p,p)\F{i}(p);
    %[Y,I] = sort(p);
    %unew = x(I);
    
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

error_inf_matlab(i) = max(abs( error ));
res_matlab(i) = norm(Eplot{i}*(F{i} - Mc{i}*unew),inf);
u_matlab{i} = unew;

t_matlab = toc

end

w = 1;
MAX = 50;
err_inf = zeros(n_level-1,MAX);
res1 = zeros(n_level-1, MAX);
res2 = zeros(n_level-1, MAX);
err_matlab = zeros(n_level-1, MAX);
u_multigrid = cell(n_level-1,1);
u_mg_debug = cell(n_level-1,1);
for start = 1:1:n_level-1
    tic;
    [thg rg] = cart2pol(a_xcp{start},a_ycp{start});
    %V{start} = uinitialfn(thg);
    %V{start} = uexactfn(thg);
    V{start} = zeros(size(F{start}));
    %V{start} = ones(size(F{start}));
    %V{start} = rand(size(F{start})) - 0.5;
    for i = start+1:1:n_level
        V{i} = zeros(size(F{i}));
    end

    [umg_cell, err_inf(start,:), res1(start,:), res2(start,:), err_matlab(start,:), cnt] = ...
        gmg_debug(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, R, n1, n2, start, w, uexact, u_matlab, Eplot, ind, uexact_pt, MAX);

    
    u_multigrid{start} = umg_cell{cnt};
    u_mg_debug{start} = umg_cell;
    
    t_mg = toc
end


% rep_res_matlab = repmat(res_matlab, 1, 2);
% xx = [0 7];
% semilogy(xx,rep_res_matlab(1,:),'b',xx,rep_res_matlab(2,:),'r',xx,rep_res_matlab(3,:),'c', ...
%          xx,rep_res_matlab(4,:),'k',xx,rep_res_matlab(5,:),'g',xx,rep_res_matlab(6,:),'m', ...
%          xx,rep_res_matlab(7,:),'--',xx,rep_res_matlab(8,:),'r--');
% hold on

err_matlab = err_matlab(end:-1:1,:);
res1 = res1(end:-1:1,:);
res2 = res2(end:-1:1,:);
error_inf_matlab = error_inf_matlab(end:-1:1);
err_inf = err_inf(end:-1:1,:);

figure
n = 1:MAX;
n = n - 1;
if n_level == 9
    semilogy(n,err_matlab(1,:),'.-',n,err_matlab(2,:),'r.-',n,err_matlab(3,:),'c^-', ...
             n,err_matlab(4,:),'k-s',n,err_matlab(5,:),'g+--',n,err_matlab(6,:),'m-d', ...
             n,err_matlab(7,:),'b*--', n,err_matlab(8,:),'r*--');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640','N=1280')
elseif n_level == 7
    semilogy(n,err_matlab(1,:),'.-',n,err_matlab(2,:),'r.-',n,err_matlab(3,:),'c^-', ...
            n,err_matlab(4,:),'k-s',n,err_matlab(5,:),'g+--',n,err_matlab(6,:),'m-d');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320')
end
title('algebraic error comparing to the solution of matlab backslash ')
xlabel('number of vcyles')
ylabel('|u\_Vcycle-u\_Backslash|_{\infty}')

figure
n = 1:MAX;
n = n - 1;
if n_level == 9
    semilogy(n,res2(1,:),'.-',n,res2(2,:),'r.-',n,res2(3,:),'c^-', ...
             n,res2(4,:),'k-s',n,res2(5,:),'g+--',n,res2(6,:),'m-d', ...
             n,res2(7,:),'b*--', n,res2(8,:),'r*--');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640','N=1280')
elseif n_level == 7
    semilogy(n,res2(1,:),'.-',n,res2(2,:),'r.-',n,res2(3,:),'c^-', ...
            n,res2(4,:),'k-s',n,res2(5,:),'g+--',n,res2(6,:),'m-d');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320')
    
end
% semilogy(n,res(1,:),'.-',n,res(2,:),'r*-');
% legend('N=20','N=10')
title('residual |f-M*u|')
xlabel('number of vcyles')
ylabel('|residual|_{\infty}')


figure
n = 1:MAX;
n = n - 1;
if n_level == 9
    semilogy(n,res1(1,:),'.-',n,res1(2,:),'r.-',n,res1(3,:),'c^-', ...
             n,res1(4,:),'k-s',n,res1(5,:),'g+--',n,res1(6,:),'m-d', ...
             n,res1(7,:),'b*--', n,res1(8,:),'r*--');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640','N=1280')
elseif n_level == 7
    semilogy(n,res1(1,:),'.-',n,res1(2,:),'r.-',n,res1(3,:),'c^-', ...
            n,res1(4,:),'k-s',n,res1(5,:),'g+--',n,res1(6,:),'m-d');
    legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320')
    
end
% semilogy(n,res(1,:),'.-',n,res(2,:),'r*-');
% legend('N=20','N=10')

title('residual |Eplot*(f-M*u)|')
xlabel('number of vcyles')
ylabel('|residual|_{\infty}')
%title(['sin(\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])
%title(['sin(\theta)+sin(',num2str(m),'\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])

% plot error of matlab and error of different number of vcycles
figure

%err_inf_matlab = cell2mat(error_inf_matlab);

n = 1:MAX;
n = n - 1;
if n_level == 9
    semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-', ...
         n,err_inf(4,:),'k-s',n,err_inf(5,:),'g+--',n,err_inf(6,:),'m-d', ...
         n,err_inf(7,:),'b.--', n,err_inf(8,:),'r*--');
legend('N=10', 'N=20', 'N=40','N=80','N=160','N=320','N=640','N=1280')
elseif n_level == 7
    semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-', ...
        n,err_inf(4,:),'k-s',n,err_inf(5,:),'g+--',n,err_inf(6,:),'m-d');
    legend('N=10','N=20','N=40','N=80','N=160','N=320');
end
% semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
% legend('N=20','N=10')
hold on

rep_err_inf_matlab = repmat(error_inf_matlab,1,2);
xx = [0 7];
if n_level == 9
    semilogy(xx,rep_err_inf_matlab(1,:),'b',xx,rep_err_inf_matlab(2,:),'r',xx,rep_err_inf_matlab(3,:),'c', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'g',xx,rep_err_inf_matlab(6,:),'m', ...
         xx,rep_err_inf_matlab(7,:),'--',xx,rep_err_inf_matlab(8,:),'r--');
elseif n_level == 7
     semilogy(xx,rep_err_inf_matlab(1,:),'b',xx,rep_err_inf_matlab(2,:),'r',xx,rep_err_inf_matlab(3,:),'c', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'g',xx,rep_err_inf_matlab(6,:),'m');
end
fs = 12;
set(gca,'Fontsize',fs)
xlabel('\fontsize{15} number of vcyles')
ylabel('\fontsize{15} |error|_{\infty}')
title('\fontsize{15} relative errors on a circle ')
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
