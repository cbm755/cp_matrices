%% test geometric multigrid method to solve -\Delta u + u  = f in a bean-shaped domain, Neumann B.C.s
epsilon = 1;
q = 6; k = 2;
uexactfn = @(x,y) x.^q.*y.^q + x.^q + y.^q + sin(k*pi*x) + sin(k*pi*y) + cos(k*pi*x) + cos(k*pi*y);
rhsfn = @(x,y) q*(q-1)*x.^(q-2).*y.^(q-2).*(x.^2+y.^2) + q*(q-1)*(x.^(q-2) + y.^(q-2)) - k^2*pi^2*(uexactfn(x,y)-x.^q-y.^q-x.^q.*y.^q) - epsilon*uexactfn(x,y);
uxfn = @(x,y) q*x.^(q-1).*y.^q + q*x.^(q-1) + k*pi*cos(k*pi*x) - k*pi*sin(k*pi*x);
uyfn = @(x,y) q*x.^q.*y.^(q-1) + q*y.^(q-1) + k*pi*cos(k*pi*y) - k*pi*sin(k*pi*y);

% the bean curve, used for computing the normals
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

sp = cscvn(pts');

sp1 = fnder(sp);
S1.type='()';
S1.subs = {1,':'};
S2.type='()';
S2.subs = {2,':'};
% derivative of parametrisation:
xp = @(t) subsref(ppval(sp1,t), S1);
yp = @(t) subsref(ppval(sp1,t), S2);

nxfn = @(t) yp(t)' ./ sqrt(xp(t)'.^2 + yp(t)'.^2); 
nyfn = @(t) -xp(t)'./ sqrt(xp(t)'.^2 + yp(t)'.^2);

x0 = -4;
x1 = 4;
y0 = -4;
y1 = 4;

%%
% 2D example on a circle
% Construct a grid in the embedding space

%dx = 0.025;
dx = 0.00625; % grid size
dx_coarsest = 0.2;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';

dy = dx;

dim = 2;  % dimension
p = 2;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

n1 = 3;
n2 = 3;

p_f2c = 1;
p_c2f = 1;

w = 0.8;

cpf = @(x,y) cpBeanInterior(x,y);

has_boundary = true;
use_ndgrid = false;
need_param = true;

disp('building grids covering the domain... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg, a_param] = ...
    build_mg_grid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary, use_ndgrid, need_param);

n_level = length(a_band);

disp('building Laplacian matrices ... ')
L = cell(n_level,1);
for i = 1:1:n_level
   L{i} = laplacian_2d_matrix(a_x1d{i}, a_y1d{i}, order, a_band{i}, a_band{i});
   L{i} = L{i} - epsilon*speye(size(L{i}));
end

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xg, a_yg, a_band, a_bdyg, p_f2c, p_c2f);

% disp('build the matrices that evaluate v at cps of boundary function on the surface')
% Ecp_Omega_S = cell(n_level,1);
% for i = 1:1:n_level
%     Ecp_Omega_S{i} = interp2_matrix(a_x1d{i},a_y1d{i},a_xcp_S{i},a_ycp_S{i},p,a_band{i});
% end
%  
% disp('build the matrices that evaluate v at cp(bdy) on course grid using values on fine grid')
% Ecp_f2c_Omega = cell(n_level-1,1);
% for i = 1:1:n_level-1
%     Ecp_f2c_Omega{i} = interp2_matrix(a_x1d{i},a_y1d{i},a_xcp{i+1}(a_bdyg{i+1}),a_ycp{i+1}(a_bdyg{i+1}),p_f2c,a_band{i});
% end
% 
% disp('build the matrices that evaluate cp for course grid of S using values on fine grid of S')
% Ecp_f2c_S = cell(n_level-1,1);
% for i = 1:1:n_level-1
%     Ecp_f2c_S{i} = interp2_matrix(a_x1d{i},a_y1d{i},a_xcp_S{i+1},a_ycp_S{i+1},p_f2c,a_band_S{i});
% end
% 
% disp('build the matrices that evaluate boundary function at cp(bdy) for course grid of $\Omega$ using values on fine grid of S')
% Ecp_f2c_Omega_S = cell(n_level-1,1);
% for i = 1:1:n_level-1
%     Ecp_f2c_Omega_S{i} = interp2_matrix(a_x1d{i},a_y1d{i},a_xcp{i+1}(a_bdyg{i+1}),a_ycp{i+1}(a_bdyg{i+1}),p_f2c,a_band_S{i});
% end

disp('setting up rhs and allocate spaces for solns')
F = cell(n_level,1);
V = cell(n_level,1);
for i = 1:1:n_level
    F{i} = rhsfn(a_xg{i},a_yg{i});
    xcp_bdy = a_xcp{i}(a_bdyg{i}); ycp_bdy = a_ycp{i}(a_bdyg{i}); param_bdy = a_param{i}(a_bdyg{i});
    F{i}(a_bdyg{i}) = ( uxfn(xcp_bdy,ycp_bdy).*nxfn(param_bdy) + uyfn(xcp_bdy,ycp_bdy).*nyfn(param_bdy) ).* a_distg{i}(a_bdyg{i}) / a_dx{i}^2;
    V{i} = zeros(size(F{i}));
end

disp('buidling matrices to deal with boundary conditions ... ')
E_out_out = cell(n_level,1);
E_out_in = cell(n_level,1); 
a_E = cell(n_level,1);
a_Ebar = cell(n_level,1);
a_Edouble = cell(n_level,1);
a_Etriple = cell(n_level,1);
for i = 1:1:n_level
    x1d = a_x1d{i}; y1d = a_y1d{i}; band = a_band{i}; dx = a_dx{i};
    I = speye(size(L{i}));
    bdy = a_bdyg{i};
    E = interp2_matrix(x1d,y1d,a_xcp{i}(bdy),a_ycp{i}(bdy),p,band);
    cpx_bar = 2*a_xcp{i}(bdy) - a_xg{i}(bdy);
    cpy_bar = 2*a_ycp{i}(bdy) - a_yg{i}(bdy);
    Ebar = interp2_matrix(x1d,y1d,cpx_bar,cpy_bar,p,band);
    cpx_double = 2*cpx_bar - a_xcp{i}(bdy);
    cpy_double = 2*cpy_bar - a_ycp{i}(bdy); 
    Edouble = interp2_matrix(x1d,y1d,cpx_double,cpy_double,p,band);
    cpx_triple = 2*cpx_double - cpx_bar;
    cpy_triple = 2*cpy_double - cpy_bar;
    Etriple = interp2_matrix(x1d,y1d,cpx_triple,cpy_triple,p,band);
    
    % quadratic interp, 2nd order
    L_bdy = (-Ebar/2 + I(bdy,:)/2)/dx^2;

    % cubic interp, 3rd order
    % L_bdy = (Edouble/6 - Ebar + E/2 + I(bdy,:)/3) / dx^2;

    % quartic interp, 4th order
    % L_bdy = (-Etriple/12 + Edouble/2 - 1.5*Ebar + 5*E/6 + I(bdy,:)/4) /dx^2;
 
    E_out_out{i} = L_bdy(:,bdy);
    E_out_in{i} = L_bdy(:,~bdy);
    L{i}(bdy,:) = L_bdy; 
    a_E{i} = E;
    a_Ebar{i} = Ebar;
    a_Edouble{i} = Edouble;
    a_Etriple{i} = Etriple;
end 

disp('pre set-up done, start to solve ...')
error_inf_matlab = zeros(n_level-1,1);
res_matlab = zeros(n_level,1);
u_matlab = cell(n_level-1,1);
uexact = cell(n_level-1,1);
for i = 1:1:n_level-1
    tic;
    
    unew = L{i} \ F{i};
        
    t_matlab = toc
    
    uexact{i} = uexactfn(a_xg{i},a_yg{i});
    error = unew - uexact{i};

    error_inf_matlab(i) = max(abs( error(~a_bdyg{i}) )) / norm(uexact{i}(~a_bdyg{i}),inf);
    
    residual = F{i} - L{i}*unew;
    res_matlab(i) = norm(residual(~a_bdyg{i}),inf) / norm(F{i}(~a_bdyg{i}));
    
    u_matlab{i} = unew;

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
%     [umg, err_inf(start,:), res(start,:)] = ...
%         gmg_test(L, a_Ebar, a_Edouble, a_Etriple, E_out_out, E_out_in, Ecp_Omega_S, Ecp_f2c_Omega, Ecp_f2c_S, Ecp_f2c_Omega_S, V, F, FonS, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, start, w, uexact, MAX);
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
%title('\fontsize{15} relative residuals in the \infty-norm')
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
%xlim([0,10])