%% Test multigrid for the coupled bulk-surface Poisson problem
R = sqrt(2);

alpha = 1;
beta = 1;

exactfnv = @(x,y) beta*exp(-x.*(x-1).*y.*(y-1));
dxfnv = @(x,y) exactfnv(x,y).*(1-2*x).*y.*(y-1);
dyfnv = @(x,y) exactfnv(x,y).*(1-2*y).*x.*(x-1);
nxfn = @(x,y) x./sqrt(x.^2+y.^2);
nyfn = @(x,y) y./sqrt(x.^2+y.^2);
rhsfnv = @(x,y) - exactfnv(x,y).*( (1-2*x).^2.*y.^2.*(y-1).^2 - 2*y.*(y-1) + ...
                                 (1-2*y).^2.*x.^2.*(x-1).^2 - 2*x.*(x-1) ) + exactfnv(x,y);
                             
exactfnu = @(x,y) 1/beta * ( alpha*exactfnv(x,y) + dxfnv(x,y).*nxfn(x,y) + dyfnv(x,y).*nyfn(x,y) );
phi = @(x,y) x.^2 + y.^2 - R^2;                             
laplace_beltrami_fun = laplace_beltrami_ls2d(phi,exactfnu);
rhsfnu = @(x,y) -laplace_beltrami_fun(x,y) + exactfnu(x,y) - alpha*exactfnv(x,y) + beta*exactfnu(x,y);

x0 = -4;
x1 = 4;
y0 = -4;
y1 = 4;

%%
% Construct a grid in the embedding space
dx = 0.00625;
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

w = 0.8;

cpf = @(x,y) cpCircleInterior(x,y,R);
cpfS = @(x,y) cpCircle(x,y,R);

disp('building grids covering the domain... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg] = ...
    build_mg_grid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, true);

disp('building grids surrounding the surface... ')
[a_band_S, a_xcp_S, a_ycp_S, a_distg_S, ~, ~, ~, ~, a_xg_S, a_yg_S] = ...
    build_mg_grid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpfS, false);

n_level = length(a_band);

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xg, a_yg, a_band, p_f2c, p_c2f);
[TMf2c_S, TMc2f_S] = helper_set_TM(a_x1d, a_y1d, a_xcp_S, a_ycp_S, a_band_S, p_f2c, p_c2f);

disp('setting up rhs and allocate spaces for solns')
F = cell(n_level,1);
G = cell(n_level,1);
V = cell(n_level,1);
U = cell(n_level,1);
for i = 1:1:n_level
    %[th_x,r_x] = cart2pol(a_xg{i},a_yg{i});
    %F{i} = rhsfn(th_x,r_x);
    F{i} = rhsfnv(a_xg{i},a_yg{i});
    F{i}(a_bdyg{i}) = 0;
    G{i} = rhsfnu(a_xcp_S{i},a_ycp_S{i});
    V{i} = zeros(size(F{i}));
    U{i} = zeros(size(G{i}));
end

disp('setting up matrices')
Lv = cell(n_level,1);
Av = cell(n_level,1);
EcpVbandU = cell(n_level,1);
Lu = cell(n_level,1);
Eu = cell(n_level,1);
Au = cell(n_level,1);
EcpUbandV = cell(n_level,1);
A = cell(n_level,1);
Eoi_v = cell(n_level,1);
Eoo_v = cell(n_level,1);
Bvu = cell(n_level,1);
Buv = cell(n_level,1);
for i = 1:1:n_level
    
    Lv{i} = laplacian_2d_matrix(a_x1d{i}, a_y1d{i}, order, a_band{i}, a_band{i});
    EcpVbandU{i} = interp2_matrix(a_x1d{i}, a_y1d{i}, a_xcp{i}(a_bdyg{i}), a_ycp{i}(a_bdyg{i}), p, a_band_S{i});   %using value of u and positions of cp v

    cpxv_bar_bdy = 2*a_xcp{i}(a_bdyg{i}) - a_xg{i}(a_bdyg{i});
    cpyv_bar_bdy = 2*a_ycp{i}(a_bdyg{i}) - a_yg{i}(a_bdyg{i});
    EbarV_bdy = interp2_matrix(a_x1d{i},a_y1d{i},cpxv_bar_bdy,cpyv_bar_bdy,2,a_band{i});
    Iv = speye(size(Lv{i}));
    ng = nnz(a_bdyg{i});
    D = spdiags(a_distg{i}(a_bdyg{i}),0,ng,ng);
    Av{i} = -Lv{i} + Iv;
    A_bdy = alpha*D*(Iv(a_bdyg{i},:) + EbarV_bdy)/2 + (-EbarV_bdy + Iv(a_bdyg{i},:))/2;
    Av{i}(a_bdyg{i},:) =  A_bdy;

    Eoo_v{i} = A_bdy(:,a_bdyg{i});
    Eoi_v{i} = A_bdy(:,~a_bdyg{i});
    
    Lu{i} = laplacian_2d_matrix(a_x1d{i}, a_y1d{i}, order, a_band_S{i}, a_band_S{i});
    Eu{i} = interp2_matrix(a_x1d{i}, a_y1d{i}, a_xcp_S{i}, a_ycp_S{i}, p, a_band_S{i});
    E1u = interp2_matrix(a_x1d{i}, a_y1d{i}, a_xcp_S{i}, a_ycp_S{i}, 1, a_band_S{i});
    Iu = speye(size(Lu{i}));
    Au{i} = -(E1u*Lu{i}-4/a_dx{i}^2*(Iu-Eu{i})) + (1+beta)*Iu;
    Lu{i} = -Lu{i} + (1+beta)*Iu;
    EcpUbandV{i} = interp2_matrix(a_x1d{i}, a_y1d{i}, a_xcp_S{i}, a_ycp_S{i}, p, a_band{i});               % interpolating value of v from band of v onto positions of cp's of u

    Bvu_tmp = sparse(length(a_band{i}),length(a_band_S{i})); 
    Bvu_tmp(a_bdyg{i},:) = -beta*D*EcpVbandU{i};
    
    Buv{i} = -alpha*EcpUbandV{i};
    A{i} = [Av{i}, Bvu_tmp; Buv{i}, Au{i}];
    
    Bvu{i} = Bvu_tmp(a_bdyg{i},:);
end

A_coarsest = A{n_level};
% E1_coarsest = interp2_matrix(a_x1d{n_level}, a_y1d{n_level}, a_xcp_S{n_level}, a_ycp_S{n_level}, 1, a_band_S{n_level});
% Mu_coarsest = E1_coarsest*Lu{n_level} - 4/a_dx{n_level}^2*(speye(size(Lu{n_level}))-Eu{n_level});

disp('pre set-up done, start to solve ...')
error_inf_matlab_u = zeros(n_level-1,1);
error_inf_matlab_v = zeros(n_level-1,1);
Eplot = cell(n_level-1,1);
uexact = cell(n_level-1,1);
vexact = cell(n_level-1,1);
u_matlab = 0;
v_matlab = 0;
for i = 1:1:n_level-1
    tic;
    
    soln = A{i} \ [F{i};G{i}];
        
    t_matlab = toc
    
    v = soln(1:length(a_band{i}));
    u = soln(length(a_band{i})+1:end);
    
    if i == 3
        u_matlab = u;
        v_matlab = v;
    end
    
    %[th,r] = cart2pol(a_xg{i},a_yg{i});
    %uexact{i} = uexactfn(th,r);
    vexact{i} = exactfnv(a_xg{i},a_yg{i});
    error_v = v - vexact{i};
    error_inf_matlab_v(i) = norm( error_v(~a_bdyg{i}), inf ) / norm(vexact{i}(~a_bdyg{i}),inf);
    
    thetas = linspace(0, 2*pi, 1000)';
    r = R*ones( size(thetas) );
    [xp, yp] = pol2cart(thetas, r);
    xp = xp(:); yp = yp(:);
    Eplot{i} = interp2_matrix( a_x1d{i}, a_y1d{i}, xp, yp, p, a_band_S{i} );
    uexact{i} = exactfnu(xp,yp);
    error_u = Eplot{i}*u - uexact{i};
    error_inf_matlab_u(i) = norm( error_u, inf ) / norm(uexact{i},inf);
end
matlab_order_u = log(error_inf_matlab_u(2:end)./error_inf_matlab_u(1:end-1))/log(2);
error_inf_matlab_u = error_inf_matlab_u(end:-1:1);
matlab_order_u = matlab_order_u(end:-1:1);

matlab_order_v = log(error_inf_matlab_v(2:end)./error_inf_matlab_v(1:end-1))/log(2);
error_inf_matlab_v = error_inf_matlab_v(end:-1:1);
matlab_order_v = matlab_order_v(end:-1:1);

MAX = 10;
err_inf_u = zeros(n_level-1,MAX);
err_inf_v = zeros(n_level-1,MAX);
for start = 1:1:n_level-1
    V{start} = zeros(size(F{start}));
    for i = start+1:1:n_level
        U{i} = zeros(size(G{i}));
        V{i} = zeros(size(F{i}));
    end
    [umg, vmg, err_inf_u(start,:), err_inf_v(start,:)] = ...
        gmg(Au, Lu, Eu, A_coarsest, Av, Eoo_v, Eoi_v, Bvu, Buv, U, V, G, F, TMf2c, TMc2f, TMf2c_S, TMc2f_S, ...
            a_band, a_band_S, a_bdyg, n1, n2, start, w, Eplot, uexact, vexact, MAX);
end

err_inf_u = err_inf_u(end:-1:1,:);
err_inf_v = err_inf_v(end:-1:1,:);

figure(1)
n = 1:MAX;
n = n - 1;
if n_level == 8
    semilogy(n,err_inf_u(1,:),'o--',n,err_inf_u(2,:),'r*--',n,err_inf_u(3,:),'g+--', ...
         n,err_inf_u(4,:),'k-s',n,err_inf_u(5,:),'c^-',n,err_inf_u(6,:),'m-d', ...
            n,err_inf_u(7,:),'bx-');
    legend('N=10','N=20','N=40','N=80','N=160','N=320','N=640')
elseif n_level == 6
    semilogy(n,err_inf_u(1,:),'o--',n,err_inf_u(2,:),'r*--',n,err_inf_u(3,:),'g+--', ...
         n,err_inf_u(4,:),'k-s',n,err_inf_u(5,:),'c^-');
    legend('N=10','N=20','N=40','N=80','N=160')
elseif n_level == 4
    semilogy(n,err_inf_u(1,:),'o--',n,err_inf_u(2,:),'r*--',n,err_inf_u(3,:),'g+--');
    legend('N=10','N=20','N=40')
end
hold on
%err_inf_matlab = cell2mat(error_inf_matlab);
rep_err_inf_matlab_u = repmat(error_inf_matlab_u,1,2);
xx = [0 MAX];
if n_level == 8
    semilogy(xx,rep_err_inf_matlab_u(1,:),'b--',xx,rep_err_inf_matlab_u(2,:),'r--',xx,rep_err_inf_matlab_u(3,:),'g', ...
         xx,rep_err_inf_matlab_u(4,:),'k',xx,rep_err_inf_matlab_u(5,:),'c', xx,rep_err_inf_matlab_u(6,:),'m-', ...
            xx,rep_err_inf_matlab_u(7,:),'b-');
elseif n_level == 6
    semilogy(xx,rep_err_inf_matlab_u(1,:),'b--',xx,rep_err_inf_matlab_u(2,:),'r--',xx,rep_err_inf_matlab_u(3,:),'g', ...
         xx,rep_err_inf_matlab_u(4,:),'k',xx,rep_err_inf_matlab_u(5,:),'c');
elseif n_level == 4
     semilogy(xx,rep_err_inf_matlab_u(1,:),'b--',xx,rep_err_inf_matlab_u(2,:),'r--',xx,rep_err_inf_matlab_u(3,:),'g');
end

% semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-');
% legend('N=20','N=10')

fs = 12;
set(gca,'Fontsize',fs)
xlabel('\fontsize{15} number of v-cycles')
ylabel('\fontsize{15} ||u^h-u||_{\infty}/||u||_{\infty}')
%xlim([0,10])

figure(2)

n = 1:MAX;
n = n - 1;
if n_level == 8
    semilogy(n,err_inf_v(1,:),'o--',n,err_inf_v(2,:),'r*--',n,err_inf_v(3,:),'g+--', ...
         n,err_inf_v(4,:),'k-s',n,err_inf_v(5,:),'c^-',n,err_inf_v(6,:),'m-d', ...
            n,err_inf_v(7,:),'bx-');
    legend('N=10','N=20','N=40','N=80','N=160','N=320','N=640')
elseif n_level == 6
    semilogy(n,err_inf_v(1,:),'o--',n,err_inf_v(2,:),'r*--',n,err_inf_v(3,:),'g+--', ...
         n,err_inf_v(4,:),'k-s',n,err_inf_v(5,:),'c^-');
    legend('N=10','N=20','N=40','N=80','N=160')
elseif n_level == 4
    semilogy(n,err_inf_v(1,:),'o--',n,err_inf_v(2,:),'r*--',n,err_inf_v(3,:),'g+--');
    legend('N=10','N=20','N=40')
end
hold on
rep_err_inf_matlab_v = repmat(error_inf_matlab_v,1,2);
xx = [0 MAX];
if n_level == 8
    semilogy(xx,rep_err_inf_matlab_v(1,:),'b--',xx,rep_err_inf_matlab_v(2,:),'r--',xx,rep_err_inf_matlab_v(3,:),'g', ...
         xx,rep_err_inf_matlab_v(4,:),'k',xx,rep_err_inf_matlab_v(5,:),'c', xx,rep_err_inf_matlab_v(6,:),'m-', ...
            xx,rep_err_inf_matlab_v(7,:),'b-');
elseif n_level == 6
    semilogy(xx,rep_err_inf_matlab_v(1,:),'b--',xx,rep_err_inf_matlab_v(2,:),'r--',xx,rep_err_inf_matlab_v(3,:),'g', ...
         xx,rep_err_inf_matlab_v(4,:),'k',xx,rep_err_inf_matlab_v(5,:),'c');
elseif n_level == 4
     semilogy(xx,rep_err_inf_matlab_v(1,:),'b--',xx,rep_err_inf_matlab_v(2,:),'r--',xx,rep_err_inf_matlab_v(3,:),'g');
end

fs = 12;
set(gca,'Fontsize',fs)
xlabel('\fontsize{15} number of v-cycles')
ylabel('\fontsize{15} ||v^h-v||_{\infty}/||v||_{\infty}')
%xlim([0,10])

i = 3;
h3=figure(3);
plot2d_compdomain(v_matlab(~a_bdyg{i}), a_xg{i}(~a_bdyg{i}), a_yg{i}(~a_bdyg{i}), a_dx{i}, a_dx{i}, 3)
hold on
plot(xp',yp','-k','LineWidth',2);
axis equal;  axis tight
xlabel('\fontsize{15}x'); ylabel('\fontsize{15}y');
set(gca,'Fontsize',12)
print(h3,'-dpng','-r450','poisson-disk-v.png');

h4=figure(4);
plot2d_compdomain(u_matlab, a_xg_S{i}, a_yg_S{i}, a_dx{i}, a_dx{i}, 4)
hold on
plot(xp',yp','-k','LineWidth',2);
axis equal;  axis tight
xlabel('\fontsize{15}x'); ylabel('\fontsize{15}y');
set(gca,'Fontsize',12)
print(h4,'-dpng','-r450','poisson-disk-u.png');