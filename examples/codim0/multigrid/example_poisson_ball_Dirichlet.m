%% test geometric multigrid method to solve -\Delta u  = f in a ball, Dirichlet B.C.s
% uexactfn = @(x,y,z) exp(sin(x.*y.*z));
% rhsfn = @(x,y,z) uexactfn(x,y,z) .* ( (y.*z).^2+(z.*x).^2+(x.*y).^2 ) .* (cos(x.*y.*z).^2 - sin(x.*y.*z)); 
k = 4;
uexactfn = @(x,y,z) sin(x).*sin(y).*sin(z) + sin(k*x).*sin(k*y).*sin(k*z);
rhsfn = @(x,y,z) -3*sin(x).*sin(y).*sin(z) - 3*k^2*sin(k*x).*sin(k*y).*sin(k*z);
g_Dirichlet = @(x,y,z) uexactfn(x,y,z);


x0 = -4;
x1 = 4;
y0 = -4;
y1 = 4;
z0 = -4;
z1 = 4;

%%

dx = 0.1/8;
dx_coarsest = 0.4;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';
z1d_coarsest = (z0:dx_coarsest:z1)';

dy = dx;
dz = dx;

dim = 3;  % dimension
p = 2;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

n1 = 3;
n2 = 3;

p_f2c = 1;
p_c2f = 1;

w = 0.8;

cpf = @cpSphereInterior;

disp('building grids covering the domain... ')
[a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg] = ...
    build_mg_grid_3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf, true);

n_level = length(a_band);

disp('building Laplacian matrices ... ')
L = cell(n_level,1);
for i = 1:1:n_level
   L{i} = laplacian_3d_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i});
end

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM_3d(a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg, a_band, a_bdyg, p_f2c, p_c2f);

disp('setting up rhs and allocate spaces for solns')
F = cell(n_level,1);
V = cell(n_level,1);
for i = 1:1:n_level
    F{i} = rhsfn(a_xg{i},a_yg{i},a_zg{i});
    F{i}(a_bdyg{i}) = g_Dirichlet(a_xcp{i}(a_bdyg{i}),a_ycp{i}(a_bdyg{i}),a_zcp{i}(a_bdyg{i}));
    V{i} = zeros(size(F{i}));
end

disp('buidling matrices to deal with boundary conditions ... ')
E_out_out = cell(n_level,1);
E_out_in = cell(n_level,1); 
a_Ebar = cell(n_level,1);
a_Edouble = cell(n_level,1);
a_Etriple = cell(n_level,1);
for i = 1:1:n_level
    x1d = a_x1d{i}; y1d = a_y1d{i}; z1d = a_z1d{i}; band = a_band{i};
    I = speye(size(L{i}));
    bdy = a_bdyg{i};
    cpx_bar = 2*a_xcp{i}(bdy) - a_xg{i}(bdy);
    cpy_bar = 2*a_ycp{i}(bdy) - a_yg{i}(bdy);
    cpz_bar = 2*a_zcp{i}(bdy) - a_zg{i}(bdy);
    Ebar = interp3_matrix(x1d,y1d,z1d,cpx_bar,cpy_bar,cpz_bar,p,band);
%     cpx_double = 2*cpx_bar - a_xcp{i}(bdy);
%     cpy_double = 2*cpy_bar - a_ycp{i}(bdy);
%     cpz_double = 2*cpz_bar - a_zcp{i}(bdy); 
%     Edouble = interp3_matrix(x1d,y1d,z1d,cpx_double,cpy_double,cpz_double,p,band);
%     cpx_triple = 2*cpx_double - cpx_bar;
%     cpy_triple = 2*cpy_double - cpy_bar;
%     cpz_triple = 2*cpz_double - cpz_bar;
%     Etriple = interp3_matrix(x1d,y1d,z1d,cpx_triple,cpy_triple,cpz_triple,p,band);
    L_bdy = (I(bdy,:) + Ebar)/2;
    %L_bdy = (I(bdy,:) + 3*Ebar - Edouble) / 3;
    %L_bdy = (I(bdy,:) + 6*Ebar - 4*Edouble + Etriple) / 4;
    E_out_out{i} = L_bdy(:,bdy);
    E_out_in{i} = L_bdy(:,~bdy);
    L{i}(bdy,:) = L_bdy; 
    a_Ebar{i} = Ebar;
    a_Edouble{i} = Edouble;
%    a_Etriple{i} = Etriple;
end 

disp('pre set-up done, start to solve ...')
error_inf_matlab = zeros(n_level-1,1);
res_matlab = zeros(n_level,1);
u_matlab = cell(n_level-1,1);
uexact = cell(n_level-1,1);
for i = 1:1:n_level-1
    uexact{i} = uexactfn(a_xg{i},a_yg{i},a_zg{i});
    
%     tic;
%     
%     unew = L{i} \ F{i};
%         
%     t_matlab = toc
%     
%     error = unew - uexact{i};
% 
%     error_inf_matlab(i) = max(abs( error(~a_bdyg{i}) )) / norm(uexact{i}(~a_bdyg{i}),inf);
%     
%     residual = F{i} - L{i}*unew;
%     res_matlab(i) = norm(residual(~a_bdyg{i}),inf) / norm(F{i}(~a_bdyg{i}));
%     
%     u_matlab{i} = unew;

end
matlab_order = log(error_inf_matlab(2:end)./error_inf_matlab(1:end-1))/log(2);

error_inf_matlab = error_inf_matlab(end:-1:1);
matlab_order = matlab_order(end:-1:1);

MAX = 10;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
u_multigrid = cell(n_level-1,1);
for start = 1:1:n_level-1
    V{start} = zeros(size(F{start}));
    %V{start} = ones(size(F{start}));
    %V{start} = rand(size(F{start})) - 0.5;
    for i = start+1:1:n_level
        V{i} = zeros(size(F{i}));
    end
    [umg, err_inf(start,:), res(start,:)] = ...
         gmg(L, E_out_out, E_out_in, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, start, w, uexact, MAX);
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
    legend('N=5', 'N=10','N=20','N=40','N=80','N=160','N=320')
elseif n_level == 6
    semilogy(n,res(1,:),'o--',n,res(2,:),'r*--',n,res(3,:),'g+--', ...
             n,res(4,:),'k-s',n,res(5,:),'c^-');
    legend('N=5', 'N=10','N=20','N=40','N=80')
elseif n_level == 5
    semilogy(n,res(1,:),'o--',n,res(2,:),'r*--',n,res(3,:),'g+--', ...
        n,res(4,:),'k-s');
    legend('N=5', 'N=10','N=20','N=40') 
elseif n_level == 4
    semilogy(n,res(1,:),'o--',n,res(2,:),'r*--',n,res(3,:),'g+--');
    legend('N=5', 'N=10','N=20')    
end
% semilogy(n,res(1,:),'.-',n,res(2,:),'r*-');
% legend('N=20','N=10')
fs = 12;
set(gca,'Fontsize',fs)
%title('\fontsize{15} relative residuals in the \infty-norm')
xlabel('\fontsize{15} number of v-cycles')
ylabel('\fontsize{15} ||f^h-A^hu^h||_{\infty}/||f^h||_{\infty}')
%title('\fontsize{15} residual |Eplot*(f-A*u)|')
%xlabel('\fontsize{15} number of v-cycles')
%ylabel('\fontsize{15} |residual|_{\infty}')
%title(['sin(\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])
%title(['sin(\theta)+sin(',num2str(m),'\theta) with p=', num2str(p), ',  res = E*(f-L*v)'])
xlim([0,10])
% plot error of matlab and error of different number of vcycles
% plot error of matlab and error of different number of vcycles
figure(2)

n = 1:MAX;
n = n - 1;
if n_level == 8
    semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', ...
         n,err_inf(4,:),'k-s',n,err_inf(5,:),'c^-',n,err_inf(6,:),'m-d', ...
            n,err_inf(7,:),'bx-');
    legend('N=5', 'N=10','N=20','N=40','N=80','N=160','N=320')
elseif n_level == 6
    semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', ...
         n,err_inf(4,:),'k-s',n,err_inf(5,:),'c^-');
    legend('N=5', 'N=10','N=20','N=40','N=80')
elseif n_level == 5
    semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--', ...
         n,err_inf(4,:),'k-s');
    legend('N=5', 'N=10','N=20','N=40')
elseif n_level == 4
    semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--');
    legend('N=5', 'N=10','N=20')
end
hold on
%err_inf_matlab = cell2mat(error_inf_matlab);
rep_err_inf_matlab = repmat(error_inf_matlab,1,2);
xx = [0 MAX];
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