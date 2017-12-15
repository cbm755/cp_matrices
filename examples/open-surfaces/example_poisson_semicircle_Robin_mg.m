%% test geometric multigrid method to solve poisson equation on a semicircle
%% Exact solutions and right hand side
shift = 0;


% k = 7.5;
% rhsfn = @(th,r)  k^2*exp( cos(k*th) ).*( sin(k*th).^2 - cos(k*th) )./r.^2;
% uexactfn = @(th) exp( cos(k*th) );
% gammafn = @(th) 2+cos(th);
% % need to pay some attention to sign because of the normal derivative
% h_robin = @(th,r) uexactfn(th) + k*(-ones(size(th))).^(th/pi).*gammafn(th).*exp(cos(k*th)).*sin(k*th)./r;

k = 10;
uexactfn = @(th) cos(th)+sin(k*th);
rhsfn = @(th,r) -cos(th)./r.^2-k^2*sin(k*th)./r.^2 - shift*uexactfn(th);
gammafn = @(th) 2+cos(th);
h_robin = @(th,r) uexactfn(th) + (-ones(size(th))).^(th/pi).*gammafn(th).*sin(th)./r - k*(-ones(size(th))).^(th/pi).*gammafn(th).*cos(k*th)./r;

x0 = -3;
x1 = 3;
y0 = -3;
y1 = 3;

%%
% 2D example on a semi-circle
% Construct a grid in the embedding space

dx = 0.003125; % grid size
dx_coarsest = 0.1;   % coarsest grid size
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

radius = sqrt(3);
cpf = @(x,y) cpSemicircle(x,y,radius);  paramf = @paramSemicircle;  

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
   
   Lc{i} = Lc{i} - shift*speye(size(Lc{i}));
   Mc{i} = Mc{i} - shift*speye(size(Lc{i}));
end

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_band, a_bdyg, p_f2c, p_c2f);

disp('buidling matrices to deal with boundary conditions ... ')
E_out_out = cell(n_level,1);
E_out_in = cell(n_level,1); 
a_Ebar = cell(n_level,1);
a_Edouble = cell(n_level,1);
a_Etriple = cell(n_level,1);
a_dist = cell(n_level,1);
for i = 1:1:n_level
    x1d = a_x1d{i}; y1d = a_y1d{i}; band = a_band{i}; dx = a_dx{i};
    I = speye(size(Lc{i}));
    bdy = logical(a_bdyg{i});
    xcp_bdy = a_xcp{i}(bdy);
    ycp_bdy = a_ycp{i}(bdy);
    E = interp2_matrix(x1d,y1d,xcp_bdy,ycp_bdy,p,band);    
    % Following expressions of nx and ny are accidently true for the semicircle case. 
    nx = ones(size(xcp_bdy));
    ny = zeros(size(xcp_bdy));
    [xg,yg,a_dist{i}] = get_proj_on_conormal_2d(a_xg{i}(bdy),a_yg{i}(bdy),xcp_bdy,ycp_bdy,nx,ny);
    xg_bar = 2*xg - a_xg{i}(bdy);
    yg_bar = 2*yg - a_yg{i}(bdy);
%     xg_bar = 2*xcp_bdy - a_xg{i}(bdy);
%     yg_bar = 2*ycp_bdy - a_yg{i}(bdy);
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
    
    ng = nnz(bdy);
    D = spdiags(a_dist{i},0,ng,ng);
    [theta_bdy,_] = cart2pol(xcp_bdy,ycp_bdy);
    GAMMA = spdiags(gammafn(theta_bdy),0,ng,ng);
    % 2nd order
    M_bdy = D*(I(bdy,:) + Ebar)/2/dx^2 + GAMMA*(-Ebar/2 + I(bdy,:)/2)/dx^2;

    % 3rd order
    % M_bdy = D*(I(bdy,:) + 3*Ebar - Edouble) / (3*dx^2) + GAMMA*(Edouble/6 - Ebar + E/2 + I(bdy,:)/3) / dx^2;
    
    E_out_out{i} = M_bdy(:,bdy);
    E_out_in{i} = M_bdy(:,~bdy);
    Mc{i}(bdy,:) = M_bdy; 
    a_Ebar{i} = Ebar;
    a_Edouble{i} = Edouble;
    a_Etriple{i} = Etriple;

end 

disp('setting up rhs and allocate spaces for solns')
F = cell(n_level,1);
V = cell(n_level,1);
for i = 1:1:n_level
    [th, r] = cart2pol(a_xcp{i},a_ycp{i});
    F{i} = rhsfn(th,r);
    bdyg = logical(a_bdyg{i});
    F{i}(bdyg) = h_robin(th(bdyg), r(bdyg)) .* a_dist{i}  / a_dx{i}^2;
    V{i} = zeros(size(F{i}));
end


disp('set up sample points and interp matrices to evaluate the errors')
% plotting grid on a semi-circle, using theta as a parameterization
thetas = linspace(0, pi, 1000)';
r = radius*ones( size(thetas) );
uexact = uexactfn(thetas);
% plotting grid in Cartesian coords
[xp, yp] = pol2cart(thetas, r);
xp = xp(:); yp = yp(:);

Eplot = cell(n_level-1,1);
for i = 1:1:n_level-1
    Eplot{i} = interp2_matrix( a_x1d{i}, a_y1d{i}, xp, yp, p, a_band{i} );
end

disp('pre set-up done, start to solve ...')
error_inf_matlab = zeros(n_level-1,1);
res_matlab = zeros(n_level,1);
u_matlab = cell(n_level-1,1);
for i = 1:1:n_level-1
    tic;
    
    unew = Mc{i} \ F{i};
        
    t_matlab = toc
    
    th = cart2pol(a_xcp{i},a_ycp{i});

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
    legend('N=20','N=40','N=80','N=160','N=320')
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
    legend('N=20','N=40','N=80','N=160','N=320')
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
%xlim([0,10])