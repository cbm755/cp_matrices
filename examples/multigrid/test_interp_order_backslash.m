
function [err_inf_matlab] = test_interp_order_backslash(p)

disp(['p = ', num2str(p)])

x0 = -3;
x1 = 3;
y0 = -3;
y1 = 3;

dx = 0.003125; % grid size
dx_coarsest = 0.2;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';


dim = 2;  % dimension
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

n1 = 2;
n2 = 1;

p_mg = 3;

w = 1;

cpf = @cpCircle;

%  rhsfn = @(th) exp( cos(th) ).*( sin(th).^2 - cos(th) );
%  uexactfn = @(th) exp( cos(th) ) - 2*exp(1);
%  rhsfn = @(th) -sin(th) - 4*sin(2*th) - 9*sin(3*th) - 16*sin(4*th) - 25*sin(5*th);
%  uexactfn = @(th) sin(th) + sin(2*th) + sin(3*th) + sin(4*th) + sin(5*th);
%  rhsfn = @(th) -cos(th);
%  uexactfn = @(th) cos(th)-2;
%  rhsfn = @(th) -sin(th);
%  uexactfn = @(th) sin(th);
%  rhsfn = @(th) -81*sin(9*th);
%  uexactfn = @(th) sin(9*th);
rhsfn = @(th) -sin(th) - 121*sin(11*th);
uexactfn = @(th) sin(th) + sin(11*th);

has_boundary = false;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%    helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);

disp('building cp grids ... ')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d] = ...
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);

disp('building cp matrices ... ')
[Mc, Lc, Ec] = build_mg_cpmatrix(a_band, a_xcp, a_ycp, a_x1d, a_y1d, p, order);

disp('building right hand side and allocate space for solution ... ')
[V, F] = helper_set_rhs(a_xcp, a_ycp, rhsfn);

% disp('building transform matrices to do restriction and prolongation later ... ')
% [TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_band, a_bdyg, p, p_mg);

disp('making the problem to solve poisson equation on a circle well-posedness ... ')
% Following line of code seems a bit confusing because a circle does not have a
% boundary, and the logical variable 'has_boundary' is set to be 'false',
% but in order to make the problem wellposedness, we do something with the
% coefficient matrix just as when dealing with Neumann Boundary Conditions,
% so we add the following line:
[Mc, Lc, Ec] = app_bnd(Mc, Lc, Ec, a_xcp, a_ycp, a_bdyg, 'neumann');



cpxg = a_xcp{1};
cpyg = a_ycp{1};
[thg, rg] = cart2pol(cpxg,cpyg);
ue = uexactfn(thg);

n_level = length(a_band);
Eplot = cell(n_level-1,1);

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0, 2*pi, 1000);
r = ones( size(thetas) );
% plotting grid in Cartesian coords
[xp, yp] = pol2cart(thetas, r);
xp = xp(:); yp = yp(:);
uexact = uexactfn(thetas);
    
dx_tmp = dx;

for i = 1:1:n_level-1
    x = (x0:dx_tmp:x1)';
    y = (y0:dx_tmp:y1)';
    
    Eplot{i} = interp2_matrix_band( x, y, xp, yp, p, a_band{i} );
    dx_tmp = 2*dx_tmp;
end

semicircplot = cell(n_level-1,1);
error_inf_matlab = cell(n_level-1,1);
for i = 1:1:n_level-1
unew = Mc{i} \ F{i};
semicircplot{i} = Eplot{i}*unew;
error_inf_matlab{i} = max(abs( uexactfn(thetas) - semicircplot{i}' ));
end

% MAX = 20;
% err_inf = zeros(n_level-1,MAX);
% for start = 1:1:n_level-1
%     [umg err_inf(start,:)] = ...
%         gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, start, w, uexact, Eplot, MAX);
% end

err_inf_matlab = cell2mat(error_inf_matlab);
err_inf_matlab = err_inf_matlab';

end

