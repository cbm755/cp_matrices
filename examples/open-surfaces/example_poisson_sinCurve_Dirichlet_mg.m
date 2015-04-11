%% Poisson equation on a parametrized sin curve with Dirichlet B.C.s
%% Exact solutions and right hand side

xs = @(x) x;
ys = @(x) sin(x);

xp = @(x) ones(size(x));
yp = @(x) cos(x);

xpp = @(x) zeros(size(x));
ypp = @(x) -sin(x);

endpt = [0, 3*pi/2];

x0 = -1;
x1 = 6;
y0 = -1.5;
y1 = 1.5;

%%
% 2D example on a semi-circle
% Construct a grid in the embedding space

dx = 0.003125/4; % grid size
dx_coarsest = 0.025;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';

dy = dx;

x1d = (x0:dx:x1)';
y1d = (y0:dx:y1)';

dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

cpf = @(x,y) cpParamCurveOpen_tmp(x,y,xs,ys,xp,yp,xpp,ypp,endpt);  
% nxfn = @(x,y) -cos(x)./sqrt(1+cos(x).^2); 
% nyfn = @(x,y) 1./sqrt(1+cos(x).^2);

% If the curve or surface has 'real' boundary:
has_boundary = true;

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg] = ...
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);

n_level = length(a_band);

disp('building Laplacian matrices ... ')
Lc = cell(n_level,1);
Mc = cell(n_level,1);
Ec = cell(n_level,1);
for i = 1:1:n_level
   Lc{i} = laplacian_2d_matrix(a_x1d{i}, a_y1d{i}, order, a_band{i}, a_band{i});
   E1 = interp2_matrix(a_x1d{i},a_y1d{i},a_xcp{i},a_ycp{i},1,a_band{i});
   E3 = interp2_matrix(a_x1d{i},a_y1d{i},a_xcp{i},a_ycp{i},3,a_band{i});
   Mc{i} = E1*Lc{i} - 2*dim/a_dx{i}^2*(speye(size(Lc{i})) - E3);
   Ec{i} = E3;
end

% parameters of v-cycle
n1 = 5;
n2 = 5;

p_f2c = 1;
p_c2f = 1;

w = 1;

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_band, a_bdyg, p_f2c, p_c2f);

% exact solutions
k = 10;
uexactfn = @(x) cos(x) + sin(k*x);
rhsfn = @(x) -(cos(x)+k^2*sin(k*x))./(1+cos(x).^2) + (-sin(x)+k*cos(k*x)).*sin(x).*cos(x)./(1+cos(x).^2).^2;

disp('setting up rhs and allocate spaces for solns')
F = cell(n_level,1);
V = cell(n_level,1);
for i = 1:1:n_level
    F{i} = rhsfn(a_xcp{i});
    bdyg = logical(a_bdyg{i});
   
    F{i}(bdyg) = uexactfn(a_xcp{i}(bdyg));
    V{i} = zeros(size(F{i}));
end

disp('buidling matrices to deal with boundary conditions ... ')
E_out_out = cell(n_level,1);
E_out_in = cell(n_level,1); 
a_Ebar = cell(n_level,1);
a_Edouble = cell(n_level,1);
a_Etriple = cell(n_level,1);
for i = 1:1:n_level
    x1d = a_x1d{i}; y1d = a_y1d{i}; band = a_band{i};
    I = speye(size(Lc{i}));
    bdy = logical(a_bdyg{i});
    xg_bar = 2*a_xcp{i}(bdy) - a_xg{i}(bdy);
    yg_bar = 2*a_ycp{i}(bdy) - a_yg{i}(bdy);
    [cpx_bar,cpy_bar] = cpf(xg_bar,yg_bar);
    Ebar = interp2_matrix(x1d,y1d,cpx_bar,cpy_bar,p,band);
    xg_double = 2*xg_bar - a_xcp{i}(bdy);
    yg_double = 2*yg_bar - a_ycp{i}(bdy); 
    [cpx_double, cpy_double] = cpf(xg_double,yg_double);
    Edouble = interp2_matrix(x1d,y1d,cpx_double,cpy_double,p,band);
    xg_triple = 2*xg_double - xg_bar;
    yg_triple = 2*yg_double - yg_bar;
    [cpx_triple, cpy_triple] = cpf(xg_triple,yg_triple);
    Etriple = interp2_matrix(x1d,y1d,cpx_triple,cpy_triple,p,band);
    M_bdy = (I(bdy,:) + Ebar)/2;
    %M_bdy = (I(bdy,:) + 3*Ebar - Edouble) / 3;
    %M_bdy = (I(bdy,:) + 6*Ebar - 4*Edouble + Etriple) / 4;
    E_out_out{i} = M_bdy(:,bdy);
    E_out_in{i} = M_bdy(:,~bdy);
    Mc{i}(bdy,:) = M_bdy; 
    a_Ebar{i} = Ebar;
    a_Edouble{i} = Edouble;
    a_Etriple{i} = Etriple;
end 

disp('set up sample points and interp matrices to evaluate the errors')
% plotting grid on a sine curve, using x as a parameterization
xs = linspace(endpt(1),endpt(2),500)';
ys = sin(xs);
uexact = uexactfn(xs);
Eplot = cell(n_level-1,1);
for i = 1:1:n_level-1
    Eplot{i} = interp2_matrix( a_x1d{i}, a_y1d{i}, xs, ys, p, a_band{i} );
end

disp('pre set-up done, start to solve ...')
error_inf_matlab = zeros(n_level-1,1);
res_matlab = zeros(n_level,1);
u_matlab = cell(n_level-1,1);
for i = 1:1:n_level-1
    tic;
    
    unew = Mc{i} \ F{i};
        
    t_matlab = toc

    error_inf_matlab(i) = norm(Eplot{i}*unew-uexact,inf) / norm(uexact,inf);
 
    u_matlab{i} = unew;

end
matlab_order = log(error_inf_matlab(2:end)./error_inf_matlab(1:end-1))/log(2);


MAX = 10;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
res2 = zeros(n_level-1, MAX);
err_matlab = zeros(n_level-1, MAX);
u_multigrid = cell(n_level-1,1);
u_mg_debug = cell(n_level-1,1);

R = cell(n_level,1);

for start = 1:1:n_level-1
    V{start} = zeros(size(F{start}));
    %V{start} = ones(size(F{start}));
    %V{start} = rand(size(F{start})) - 0.5;
    for i = start+1:1:n_level
        V{i} = zeros(size(F{i}));
    end
    [umg, err_inf(start,:), res(start,:)] = ...
        gmg(Mc, Lc, Ec, E_out_out, E_out_in, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, start, w, Eplot, uexact, MAX);
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
    legend('N=80','N=160','N=320','N=640','N=1280')
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
xlim([0,10])

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
    legend('N=80','N=160','N=320','N=640','N=1280')
elseif n_level == 4
    semilogy(n,err_inf(1,:),'o--',n,err_inf(2,:),'r*--',n,err_inf(3,:),'g+--');
    legend('N=10','N=20','N=40')
end
hold on
%err_inf_matlab = cell2mat(error_inf_matlab);
error_inf_matlab = error_inf_matlab(end:-1:1);
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