%% Test multigrid for the coupled bulk-surface Poisson problem
R = 1;

alpha = 1;
beta = 1;

tic
disp('setting up exact solutions and right hand side functions')
% k = 4;
% exactfnv = @(x,y,z) sin(x).*sin(y).*sin(z) + sin(k*x).*sin(k*y).*sin(k*z);
% rhsfnv = @(x,y,z) -3*sin(x).*sin(y).*sin(z) - 3*k^2*sin(k*x).*sin(k*y).*sin(k*z);
% dxfnv = @(x,y,z) cos(x).*sin(y).*sin(z) + k*cos(k*x).*sin(k*y).*sin(k*z);
% dyfnv = @(x,y,z) sin(x).*cos(y).*sin(z) + k*sin(k*x).*cos(k*y).*sin(k*z);
% dzfnv = @(x,y,z) sin(x).*sin(y).*cos(z) + k*sin(k*x).*sin(k*y).*cos(k*z);

exactfnv = @(x,y,z) beta*exp(-x.*(x-1).*y.*(y-1).*z.*(z-1));
dxfnv = @(x,y,z) exactfnv(x,y,z).*(1-2*x).*y.*(y-1).*z.*(z-1);
dyfnv = @(x,y,z) exactfnv(x,y,z).*(1-2*y).*x.*(x-1).*z.*(z-1);
dzfnv = @(x,y,z) exactfnv(x,y,z).*(1-2*z).*x.*(x-1).*y.*(y-1);
rhsfnv = @(x,y,z) - exactfnv(x,y,z).*( ((1-2*x).*y.*(y-1).*z.*(z-1)).^2 - 2*y.*(y-1).*z.*(z-1) + ...
                                     ((1-2*y).*x.*(x-1).*z.*(z-1)).^2 - 2*x.*(x-1).*z.*(z-1) + ...
                                     ((1-2*z).*x.*(x-1).*y.*(y-1)).^2 - 2*x.*(x-1).*y.*(y-1) ) ... 
                  + exactfnv(x,y,z);
              
nxfn = @(x,y,z) x./sqrt(x.^2+y.^2+z.^2);
nyfn = @(x,y,z) y./sqrt(x.^2+y.^2+z.^2);
nzfn = @(x,y,z) z./sqrt(x.^2+y.^2+z.^2);

                             
exactfnu = @(x,y,z) 1/beta * ( alpha*exactfnv(x,y,z) + dxfnv(x,y,z).*nxfn(x,y,z) + dyfnv(x,y,z).*nyfn(x,y,z) + dzfnv(x,y,z).*nzfn(x,y,z) );
phi = @(x,y,z) x.^2 + y.^2 + z.^2 - R^2;                             
laplace_beltrami_fun = laplace_beltrami_ls3d(phi,exactfnu);
rhsfnu = @(x,y,z) -laplace_beltrami_fun(x,y,z) + exactfnu(x,y,z) - alpha*exactfnv(x,y,z) + beta*exactfnu(x,y,z);
disp('done')
toc

x0 = -4;
x1 = 4;
y0 = x0;
y1 = x1;
z0 = x0;
z1 = x1;

%%
% Construct a grid in the embedding space
dx = 0.1/2;
dx_coarsest = 0.4;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';
z1d_coarsest = (z0:dx_coarsest:z1)';

dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

n1 = 3;
n2 = 3;

p_f2c = 1;
p_c2f = 1;

w = 0.8;

cpf = @(x,y,z) cpSphereInterior(x,y,z,R);
cpfS = @(x,y,z) cpSphere(x,y,z,R);

disp('building grids covering the domain... ')
[a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg] = ...
    build_mg_grid_3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf, true);

disp('building grids surrounding the surface... ')
[a_band_S, a_xcp_S, a_ycp_S, a_zcp_S, a_distg_S, ~, ~, ~, ~, ~, a_xg_S, a_yg_S, a_zg_S] = ...
    build_mg_grid_3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpfS, false);

n_level = length(a_band);

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM_3d(a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg, a_band, p_f2c, p_c2f);
[TMf2c_S, TMc2f_S] = helper_set_TM_3d(a_x1d, a_y1d, a_z1d, a_xcp_S, a_ycp_S, a_zcp_S, a_band_S, p_f2c, p_c2f);

disp('setting up rhs and allocate spaces for solns')
F = cell(n_level,1);
G = cell(n_level,1);
V = cell(n_level,1);
U = cell(n_level,1);
for i = 1:1:n_level
    %[th_x,r_x] = cart2pol(a_xg{i},a_yg{i});
    %F{i} = rhsfn(th_x,r_x);
    F{i} = rhsfnv(a_xg{i},a_yg{i},a_zg{i});
    F{i}(a_bdyg{i}) = 0;
    G{i} = rhsfnu(a_xcp_S{i},a_ycp_S{i},a_zcp_S{i});
    V{i} = zeros(size(F{i}));
    U{i} = zeros(size(G{i}));
end
disp('done')

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
    
    Lv{i} = laplacian_3d_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i});
    %using value of u and positions of cp v
    EcpVbandU{i} = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}(a_bdyg{i}), a_ycp{i}(a_bdyg{i}), a_zcp{i}(a_bdyg{i}), p, a_band_S{i});   

    cpxv_bar_bdy = 2*a_xcp{i}(a_bdyg{i}) - a_xg{i}(a_bdyg{i});
    cpyv_bar_bdy = 2*a_ycp{i}(a_bdyg{i}) - a_yg{i}(a_bdyg{i});
    cpzv_bar_bdy = 2*a_zcp{i}(a_bdyg{i}) - a_zg{i}(a_bdyg{i});
    EbarV_bdy = interp3_matrix(a_x1d{i},a_y1d{i},a_z1d{i},cpxv_bar_bdy,cpyv_bar_bdy,cpzv_bar_bdy,2,a_band{i});
    Iv = speye(size(Lv{i}));
    ng = nnz(a_bdyg{i});
    D = spdiags(a_distg{i}(a_bdyg{i}),0,ng,ng);
    Av{i} = -Lv{i} + Iv;
    A_bdy = alpha*D*(Iv(a_bdyg{i},:) + EbarV_bdy)/2 + (-EbarV_bdy + Iv(a_bdyg{i},:))/2;
    Av{i}(a_bdyg{i},:) =  A_bdy;

    Eoo_v{i} = A_bdy(:,a_bdyg{i});
    Eoi_v{i} = A_bdy(:,~a_bdyg{i});
    
    Lu{i} = laplacian_3d_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band_S{i}, a_band_S{i});
    Eu{i} = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp_S{i}, a_ycp_S{i}, a_zcp_S{i}, p, a_band_S{i});
    E1u = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp_S{i}, a_ycp_S{i}, a_zcp_S{i}, 1, a_band_S{i});
    Iu = speye(size(Lu{i}));
    Au{i} = -(E1u*Lu{i}-6/a_dx{i}^2*(Iu-Eu{i})) + (1+beta)*Iu;
    Lu{i} = -Lu{i} + (1+beta)*Iu;
    % interpolating value of v from band of v onto positions of cp's of u
    EcpUbandV{i} = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp_S{i}, a_ycp_S{i}, a_zcp_S{i}, p, a_band{i});               
    
    Buv{i} = -alpha*EcpUbandV{i};
    
    if abs( i - n_level ) <= 2
        Bvu_tmp = sparse(length(a_band{i}),length(a_band_S{i})); 
        Bvu_tmp(a_bdyg{i},:) = -beta*D*EcpVbandU{i};
        A{i} = [Av{i}, Bvu_tmp; Buv{i}, Au{i}];
    end
    
    Bvu{i} = -beta*D*EcpVbandU{i};
end
disp('done')

A_coarsest = A{n_level};
% E1_coarsest = interp2_matrix(a_x1d{n_level}, a_y1d{n_level}, a_xcp_S{n_level}, a_ycp_S{n_level}, 1, a_band_S{n_level});
% Mu_coarsest = E1_coarsest*Lu{n_level} - 4/a_dx{n_level}^2*(speye(size(Lu{n_level}))-Eu{n_level});


error_inf_matlab_u = zeros(n_level-1,1);
error_inf_matlab_v = zeros(n_level-1,1);
Eplot = cell(n_level-1,1);
uexact = cell(n_level-1,1);
vexact = cell(n_level-1,1);
[xp, yp, zp] = paramSphere(256,R);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
for i = 1:1:n_level-1
    vexact{i} = exactfnv(a_xg{i},a_yg{i},a_zg{i});
    Eplot{i} = interp3_matrix( a_x1d{i}, a_y1d{i}, a_z1d{i}, xp1, yp1, zp1, p, a_band_S{i} );
    uexact{i} = exactfnu(xp1,yp1,zp1);
    
    if abs(i-n_level)<=0
    tic;
    
    soln = A{i} \ [F{i};G{i}];
        
    t_matlab = toc
    
    v = soln(1:length(a_band{i}));
    u = soln(length(a_band{i})+1:end);
    
    error_v = v - vexact{i};
    error_inf_matlab_v(i) = norm( error_v(~a_bdyg{i}), inf ) / norm(vexact{i}(~a_bdyg{i}),inf);
    
    error_u = Eplot{i}*u - uexact{i};
    error_inf_matlab_u(i) = norm( error_u, inf ) / norm(uexact{i},inf);
    end
end
matlab_order_u = log(error_inf_matlab_u(2:end)./error_inf_matlab_u(1:end-1))/log(2);
error_inf_matlab_u = error_inf_matlab_u(end:-1:1);
matlab_order_u = matlab_order_u(end:-1:1);

matlab_order_v = log(error_inf_matlab_v(2:end)./error_inf_matlab_v(1:end-1))/log(2);
error_inf_matlab_v = error_inf_matlab_v(end:-1:1);
matlab_order_v = matlab_order_v(end:-1:1);

disp('pre set-up done, start to solve ...')
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
if n_level == 6
    semilogy(n,err_inf_u(1,:),'o--',n,err_inf_u(2,:),'r*--',n,err_inf_u(3,:),'g+--', ...
         n,err_inf_u(4,:),'k-s',n,err_inf_u(5,:),'c^-');
    legend('N=5','N=10','N=20','N=40','N=80')
elseif n_level  ==5
    semilogy(n,err_inf_u(1,:),'o--',n,err_inf_u(2,:),'r*--',n,err_inf_u(3,:),'g+--', ...
         n,err_inf_u(4,:),'k-s');
    legend('N=5','N=10','N=20','N=40')
elseif n_level == 4
    semilogy(n,err_inf_u(1,:),'o--',n,err_inf_u(2,:),'r*--',n,err_inf_u(3,:),'g+--');
    legend('N=5','N=10','N=20')
end
hold on
%err_inf_matlab = cell2mat(error_inf_matlab);
rep_err_inf_matlab_u = repmat(error_inf_matlab_u,1,2);
xx = [0 MAX];
if n_level == 6
    semilogy(xx,rep_err_inf_matlab_u(1,:),'b--',xx,rep_err_inf_matlab_u(2,:),'r--',xx,rep_err_inf_matlab_u(3,:),'g', ...
         xx,rep_err_inf_matlab_u(4,:),'k',xx,rep_err_inf_matlab_u(5,:),'c');
elseif n_level == 5
    semilogy(xx,rep_err_inf_matlab_u(1,:),'b--',xx,rep_err_inf_matlab_u(2,:),'r--',xx,rep_err_inf_matlab_u(3,:),'g', ...
         xx,rep_err_inf_matlab_u(4,:),'k');
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
if n_level == 6
    semilogy(n,err_inf_v(1,:),'o--',n,err_inf_v(2,:),'r*--',n,err_inf_v(3,:),'g+--', ...
         n,err_inf_v(4,:),'k-s',n,err_inf_v(5,:),'c^-');
    legend('N=5','N=10','N=20','N=40','N=80')
elseif n_level == 5
    semilogy(n,err_inf_v(1,:),'o--',n,err_inf_v(2,:),'r*--',n,err_inf_v(3,:),'g+--', ...
         n,err_inf_v(4,:),'k-s');
    legend('N=5','N=10','N=20','N=40')
elseif n_level == 4
    semilogy(n,err_inf_v(1,:),'o--',n,err_inf_v(2,:),'r*--',n,err_inf_v(3,:),'g+--');
    legend('N=5','N=10','N=20')
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

