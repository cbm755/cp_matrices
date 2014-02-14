%% test geometric multigrid method to solve poisson equation on a circle

%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');

x0 = -3;
x1 = 3;
y0 = -3;
y1 = 3;

%%
% 2D example on a circle
% Construct a grid in the embedding space

%dx = 0.003125/4; % grid size
dx = 0.003125;
%dx = 0.025;
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

n1 = 2;
n2 = 2;

p_f2c = 1;
p_c2f = 1;

w = 1;

cpf = @cpCircle;
%cpf = @cpEggCurve;

has_boundary = false;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg] = ...
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);
disp('done')

n_level = length(a_band);

disp('building cp matrices on each level of V-cycle')
[Mc, Lc, Ec, Eic, Rc, a_iband, a_oband, a_xcp, a_ycp, a_xg, a_yg, a_innerInOuter] = ...
    build_mg_ops_and_bands(a_band, a_xcp, a_ycp, a_xg, a_yg, a_x1d, a_y1d, p, order);
disp('done')

% for k = 1:1:n_level
% 
% 	% test the approach getting rid of the same rows of the matrix
% 	Mc{k} = Eic{k}*Mc{k};
% 	X = [a_xcp{k} a_ycp{k}];
% 	% finding whether the closest points of two grid points are the same
% 	[idx, dist] = rangesearch(X,X,1e-10);
% 	sidx = cellfun('size',idx,2);
% 	% if the size is bigger than 1, that indicates at least two closest
% 	% points are the same.
% 	flag = (sidx>1);
% 	tmp = idx(flag);
% 	s = nnz(flag);
%     count = 0;
% 	for i = 1:1:s
% 		j = min(tmp{i});
%         if(flag(j)==true)
%             for m = 1:1:size(tmp{i},2)
%                 n = tmp{i}(m);
%                 if (n ~= j)
%                     %[j n]
%                     count = count + 1;
%                     Mc{k}(n,:) = rand(1,size(Mc{k},1));
%                     %Mc{k}(n,:) = 0;
%                     %Mc{k}(n,j) = -1;
%                     %Mc{k}(n,n) = 1;
%                     Eic{k}(n,:) = rand(1,size(Mc{k},1));
%                     %Eic{k}(n,:) = 0;
%                     %Eic{k}(n,j) = -1;
%                     %Eic{k}(n,n) = 1;
%                 end
%             end 
%             flag(tmp{i}) = false;  
%         end
% 	end
% 
% 	%Mc{k} = lapsharp_unordered(Lc{k}, E, Rc{k});
%     size(Mc{k})
%     rank(full(Eic{k}))
%     rank(full(Mc{k}))
%     count
% end

dxc = cell(n_level,1);
dx_tmp = dx;
for i = 1:1:n_level-1
    dxc{i} = dx_tmp;
    dx_tmp = 2*dx_tmp;
end
dxc{n_level} = dx_tmp;


for k = 1:1:n_level
    E = Ec{k};
	E(a_innerInOuter{k},:) = speye(size(Ec{k},2));
    Mc{k} = Lc{k}*E;
    
    %lambda = 2*dim/dxc{i};
    %Mc{k} = Eic{k}*(Lc{k}*Ec{k}) - lambda*(speye(size(Eic{k}))-Eic{k});
    
	%Mc{k} = lapsharp_unordered(Lc{k}, E, Rc{k});
end

disp('building transform matrices to do restriction and prolongation later ... ')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_iband, a_bdyg, p_f2c, p_c2f);
disp('done')


disp('extracting diagonal entries of L, making jacobi iteration faster...')
tic;
D = cell(n_level,1);
for i = 1:1:n_level
    [i1,j1,r1] = find(Rc{i});
    Ldiagpad = Rc{i}.*Lc{i};
    D{i} = diag(Ldiagpad(i1,j1));
end
toc;
disp('done')

%% building E_plot for purpose of plotting and debug
% plotting grid on circle, using theta as a parameterization
thetas = linspace(0, 2*pi, 1000)';
r = ones( size(thetas) );
% plotting grid in Cartesian coords
[xp, yp] = pol2cart(thetas, r);
xp = xp(:); yp = yp(:);

Eplot = cell(n_level-1,1);
dx_tmp = dx;
for i = 1:1:n_level-1
    x = (x0:dx_tmp:x1)';
    y = (y0:dx_tmp:y1)';
    
    Eplot{i} = interp2_matrix_test( x, y, xp, yp, p );
    Eplot{i} = Eplot{i}(:,a_iband{i});

    dx_tmp = 2*dx_tmp;
end



%% Setting up right hand side
%rhsfn = @(th,r) ( exp( cos(th) ).*( sin(th).^2 - cos(th) ) ) ./ (r.^2);
%uexactfn = @(th) exp( cos(th) ) - exp(1);
 
%  rhsfn = @(th,r) ( 121*exp( cos(11*th) ).*( sin(11*th).^2 - cos(11*th) ) ) ./ (r.^2);
%  uexactfn = @(th) exp( cos(11*th) ) - exp(1);
% rhsfn = @(th,r) 4*exp( cos(2*th) ).*( sin(2*th).^2 - cos(2*th) )./(r.^2);
% uexactfn = @(th) exp( cos(2*th) ) - exp(1);

% rhsfn = @(th,r) ( -( 1 + 2*cos(th) ) ./ ( 2 + cos(th) ).^2 ) ./ (r.^2);
% uexactfn = @(th) log( 2 + cos(th) ) - log(3);
 
% n = 100;
% rhsfn = @(th,r) ( n*(n-1)*sin(th).^(n-2).*cos(th).^2 - n*sin(th).^n ) ./ (r.^2);
% uexactfn = @(th) sin(th).^n;
 
% n_mode = 100;
% rhsfn = @(th) rhsfn_mix_mode(th,n_mode);
% uexactfn = @(th) uexactfn_mix_mode(th,n_mode);

n = 3;
rhsfn = @(th,r) rhsfn_handle(th,r,n);
uexactfn = @(th) abs(sin(th)).^n;


%  m = 20;
%  rhsfn = @(th,r) (-sin(th) - m*sin(m*th))./(r.^2);
%  uexactfn = @(th) sin(th) + sin(m*th)/m;

%  rhsfn = @(th,r) ( -sin(th) - 4*sin(2*th) - 9*sin(3*th) - 16*sin(4*th) - 25*sin(5*th) ) ./ (r.^2);
%  uexactfn = @(th) sin(th) + sin(2*th) + sin(3*th) + sin(4*th) + sin(5*th);

% rhsfn = @(th,r) -cos(th)./(r.^2);
% uexactfn = @(th) cos(th)-1;

% rhsfn = @(th,r) -sin(th)./(r.^2);
% uexactfn = @(th) sin(th);

uexact = uexactfn(thetas);

uexact_debug = cell(n_level,1);
for i = 1:1:n_level
    [thg, rg] = cart2pol(a_xcp{i}, a_ycp{i});
    uexact_debug{i} = uexactfn(thg);
end
% v_initial_fn = @(th) sin(th);
% uexact = zeros(size(thetas));

disp('building right hand side and allocate space for solution ... ')
[V, F] = helper_set_rhs(a_xg, a_yg, rhsfn, 1);
%[V, F] = helper_set_rhs(a_xcp, a_ycp, rhsfn, 1);
disp('done')

disp('making the problem to solve poisson equation on a circle well-posedness ... ')
% Following line of code seems a bit confusing because a circle does not have a
% boundary, and the logical variable 'has_boundary' is set to be 'false',
% but in order to make the problem wellposedness, we do something with the
% coefficient matrix just as when dealing with Neumann Boundary Conditions,
% so we add the following line:
% [Mc, Lc, Ec, F] = app_bnd_2b(Mc, Lc, Ec, F, a_xcp, a_ycp, a_bdyg, 'neumann');
 [Mc, Lc, Ec, F] = app_bnd_2b(Mc, Lc, Ec, F, a_innerInOuter, a_xg, a_yg, a_bdyg, 'neumann');
disp('done')


circplot = cell(n_level-1,1);
error_inf_matlab = zeros(n_level-1,1);
res_matlab = zeros(n_level,1);
u_matlab = cell(n_level-1,1);
disp('start solving the linear system by matlab')
for i = 1:1:n_level-1
    
    %unew = Mc{i} \ F{i};
    
     m = Mc{i};
     p = amd(m);
     x = m(p,p)\F{i}(p);
     [Y,I] = sort(p);
     unew = x(I);
    
    %unew = agmg(Mc{i},F{i},[],1e-12,1000);
    
 %unew = Ec{i}*unew;
% By increasing the maximal number of iterations, bicgstab will converge;
% however as grid become finer, number of iterations will increase sharply
% [unew flag] = bicgstab(Mc{i}, F{i}, 1e-10, 200);

% gmres seems not converge for this problem
% unew = gmres(Mc{i}, F{i}, 3, 1e-10);

circplot{i} = Eplot{i}*unew;
error_inf_matlab(i) = max(abs( uexact - circplot{i} ));
res_matlab(i) = norm(Eplot{i}*(F{i} - Mc{i}*unew),inf);

u_matlab{i} = unew;

end
disp('done')

MAX = 50;
err_inf = zeros(n_level-1,MAX);
res = zeros(n_level-1, MAX);
u_multigrid = cell(n_level-1,1);
for start = 1:1:n_level-1
   [umg err_inf(start,:) res(start,:)] = ...
       gmg_2b(Mc, Lc, Ec, Eic, D, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, start, w, uexact, Eplot, MAX);
%   [umg err_inf(start,:) res(start,:)] = gmg_M(Mc, Lc, Eic, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, start, w, uexact, Eplot, MAX); 
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

if n_level == 9
    semilogy(n,res(1,:),'.-',n,res(2,:),'r.-',n,res(3,:),'c^-', ...
             n,res(4,:),'k-s',n,res(5,:),'g+--',n,res(6,:),'m-d', ...
             n,res(7,:),'b*--', n,res(8,:),'r*--');
    legend('N=1280', 'N = 640', 'N=320','N=160','N=80','N=40','N=20','N=10')
elseif n_level == 7
    semilogy(n,res(1,:),'.-',n,res(2,:),'r.-',n,res(3,:),'c^-', ...
            n,res(4,:),'k-s',n,res(5,:),'g+--',n,res(6,:),'m-d');
    legend('N=320','N=160','N=80','N=40','N=20','N=10')
    
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
if n_level == 9
    semilogy(xx,rep_err_inf_matlab(1,:),'b',xx,rep_err_inf_matlab(2,:),'r',xx,rep_err_inf_matlab(3,:),'c', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'g',xx,rep_err_inf_matlab(6,:),'m', ...
         xx,rep_err_inf_matlab(7,:),'--',xx,rep_err_inf_matlab(8,:),'r--');
elseif n_level == 7
     semilogy(xx,rep_err_inf_matlab(1,:),'b',xx,rep_err_inf_matlab(2,:),'r',xx,rep_err_inf_matlab(3,:),'c', ...
         xx,rep_err_inf_matlab(4,:),'k',xx,rep_err_inf_matlab(5,:),'g',xx,rep_err_inf_matlab(6,:),'m');
end
hold on

n = 1:MAX;
if n_level == 9
    semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-', ...
         n,err_inf(4,:),'k-s',n,err_inf(5,:),'g+--',n,err_inf(6,:),'m-d', ...
         n,err_inf(7,:),'b.--', n,err_inf(8,:),'r*--');
legend('N=1280', 'N = 640', 'N=320','N=160','N=80','N=40','N=20','N=10')
elseif n_level == 7
    semilogy(n,err_inf(1,:),'.-',n,err_inf(2,:),'r*-',n,err_inf(3,:),'c^-', ...
        n,err_inf(4,:),'k-s',n,err_inf(5,:),'g+--',n,err_inf(6,:),'m-d');
    legend('N=320','N=160','N=80','N=40','N=20','N=10');
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

