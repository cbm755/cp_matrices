%% Heat equation on a circle
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

% This example demonstrates two bands as in the implicit CP paper
% [Macdonald, Ruuth 2009]


%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../../cp_matrices');

% add functions for finding the closest points
addpath('../../surfaces');

x0 = -2;
x1 = 2;
y0 = -2;
y1 = 2;

%%
% 2D example on a circle
% Construct a grid in the embedding space

%dx = 0.003125; % grid size
%dx = 0.0125;
dx = 0.1
dx_coarsest = 0.2;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';

dy = dx;
dt = dx/4;

dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

w = 1;
n1 = 3;
n2 = 3;
nc = 3;

MAX = 2;

p_f2c = 2;
p_c2f = 2;

cpf = @cpCircle;
uexactfn = @(t,th) exp(-t)*cos(th) + exp(-9*t)*cos(3*th) + exp(-25*t)*cos(5*th) + ...
              exp(-4*t)*sin(2*th) + exp(-16*t)*sin(4*th) + exp(-36*t)*sin(6*th);
%uexactfn = @(t,th) exp(-t)*cos(th) + exp(-9*t)*cos(3*th);
%uexactfn = @(t,th) exp(-t)*sin(th) + exp(-9*t)*sin(3*th);
%uexactfn = @(t,th) exp(-t)*cos(th);
%uexactfn = @(t,th) exp(-t)*sin(th);

has_boundary = false;

%rhsfn = @(th) uexactfn(0,th);
%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%       helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);
disp('building cp grids on each level of V-cycle')
[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d] = ...
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary, 0, 0);
disp('done')
disp('building cp matrices on each level of V-cycle')
[Mc, Lc, Ec] = build_mg_cpmatrix(a_band, a_xcp, a_ycp, a_x1d, a_y1d, p, order);
disp('done')

n_level = length(a_band);
V = cell(n_level,1);
F = cell(n_level,1);
A_Mc = cell(n_level,1);
A_Lc = cell(n_level,1);
for i = 1:1:n_level
    A_Mc{i} = speye(size(Mc{i})) - 0.5*dt*Mc{i};
    A_Lc{i} = speye(size(Lc{i})) - 0.5*dt*Lc{i}; 
    %A_Mc{i} = speye(size(Mc{i})) - dt*Mc{i};
    %A_Lc{i} = speye(size(Lc{i})) - dt*Lc{i}; 
end 


disp('building transform matrices doing restriction and prolongation ...')
[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_band, a_bdyg, p_f2c, p_c2f);
disp('done')

cpxg = a_xcp{1};
cpyg = a_ycp{1};
[thg, rg] = cart2pol(cpxg,cpyg);
%u0 = cos(thg) + cos(3*thg);
%u0 = sin(thg) + sin(3*thg);
%u0 = cos(thg);
u0 = uexactfn(0,thg);

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,1000)';
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);

Eplot = cell(n_level-1, 1);
dx_tmp = dx;
for i = 1:1:n_level-1
    x = (x0:dx_tmp:x1)';
    y = (y0:dx_tmp:y1)';

    Eplot{i} = interp2_matrix( x, y, xp, yp, p );
    Eplot{i} = Eplot{i}(:, a_band{i});
    dx_tmp = 2*dx_tmp;
end


%% Time-stepping for the heat equation
Tf = dt;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;

gap = 200;
err_plot = zeros(numtimesteps,1);
err_plot1 = zeros(numtimesteps,1);
t_plot = dt:dt:Tf;

figure(1);
u = u0;
tic
for kt = 1:numtimesteps
    
    t = kt*dt;
    %unew = Mc{1} \ u;
    
    % Both bicgstab and gmres works well in this case, and seems to be
    % a bit faster than multigrid; and I think gmres works best of the
    % three
    %[unew flag relres iter] = bicgstab(A_Mc{1}, u, 1e-10, 20, [], [], u);
    %[unew flag relres iter] = gmres(A_Mc{1}, u, 10, 1e-10, 20, [], [], u);
	u = u + 0.5*dt*Mc{1}*u;
    unew = A_Mc{1} \ u;
    %[unew flag] = gmres(A_Mc{1}, u, 5, 1e-6);
    res = Eplot{1}*(u - A_Mc{1}*unew);
    r = norm(res,inf)
    %u = Ec{1}*unew;
    u = unew;
    
    circplot = Eplot{1}*u;
    error_circ_inf = max(abs( uexactfn(t,thetas) - circplot )) / max(abs(uexactfn(t,thetas)));
    err_plot1(kt) = error_circ_inf;
    %plotting
%     if ( (kt < 5) || (mod(kt,gap) == 0) || (kt == numtimesteps) )
%       t
%       % plot in the embedded domain: shows the computational band
%       plot2d_compdomain(u, a_xcp{1}, a_ycp{1}, dx, dy, 1);
%       hold on;
%       plot(xp,yp,'k-', 'linewidth', 2);
%       title( ['embedded domain: soln at time ' num2str(t) ...
%             ', timestep #' num2str(kt)] );
% 
%       % plot value on circle
%       figure(2); clf;
%       
%       plot(thetas, circplot);
%       title( ['soln at time ' num2str(t) ', on circle'] );
%       xlabel('theta'); ylabel('u');
%       hold on;
%       % plot analytic result
%       plot(thetas, uexactfn(t,thetas), 'r--');
%       plot(thetas, Eplot{1}*u0, 'g-.');
%       legend('explicit Euler', 'exact answer', 'initial condition ', ...
%            'Location', 'SouthEast');
%       
%       [dx dt t error_circ_inf]
%        
%       pause(0);
%     end
end
time_matlab = toc

close all
u = u0;
figure(1);
tic
for kt = 1:numtimesteps

    t = kt*dt;
    %% implicit backward euler
    % set the initial value of multigrid as the solution of last time step
    %[V, F] = helper_set_V_F(a_band, u);
	V{1} = u;
	for i = 2:1:n_level
		V{i} = zeros(size(a_band{i}));
	end
	F{1} = u + 0.5*dt*Lc{1}*u;
    %F{1} = u;
    %unew = gmg_M(A_Mc, A_Lc, Ec, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w, uexactfn(t,thetas), Eplot, MAX);
    unew = gmg(A_Mc, A_Lc, Ec, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w, uexactfn(t,thetas), Eplot, MAX);
    %unew = helper_vcycle(A_Mc, A_Lc, Ec, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w);
    res = Eplot{1}*(u - A_Mc{1}*unew);
    r = norm(res,inf)
    
    u = unew;

    circplot = Eplot{1}*u;
    error_circ_inf = max(abs( uexactfn(t,thetas) - circplot )) / max(abs(uexactfn(t,thetas)));
    err_plot(kt) = error_circ_inf;
    %plotting
%     if ( (kt < 5) || (mod(kt,gap) == 0) || (kt == numtimesteps) )
%       t
%       % plot in the embedded domain: shows the computational band
%       plot2d_compdomain(u, a_xcp{1}, a_ycp{1}, dx, dx, 1);
%       hold on;
%       plot(xp,yp,'k-', 'linewidth', 2);
%       title( ['embedded domain: soln at time ' num2str(t) ...
%             ', timestep #' num2str(kt)] );
% 
%       % plot value on circle
%       figure(2); clf;
%       
%       plot(thetas, circplot);
%       title( ['soln at time ' num2str(t) ', on circle'] );
%       xlabel('theta'); ylabel('u');
%       hold on;
%       % plot analytic result
%       plot(thetas, uexactfn(t,thetas), 'r--');
%       plot(thetas, Eplot{1}*u0, 'g-.');
%       legend('explicit Euler', 'exact answer', 'initial condition ', ...
%            'Location', 'SouthEast');
%       
%       [dx dt t error_circ_inf]
%        
%       pause(0);
%     end
end
time_mg = toc

figure(3)
plot(t_plot, err_plot1, t_plot,err_plot,'r')
xlabel('time')
ylabel('relative error')
%title(['relative error of e^{-t}cos(\theta)+e^{-9t}cos(3\theta),  dx = ', num2str(dx)])
title(['relative error of e^{-t}sin(\theta)+e^{-9t}sin(3\theta),  dx = ', num2str(dx)])
legend('matlab','multigrid')



E1 = interp2_matrix(a_x1d{1}, a_y1d{1}, a_xcp{1}, a_ycp{1}, 1,a_band{1});
E3 = interp2_matrix(a_x1d{1}, a_y1d{1}, a_xcp{1}, a_ycp{1}, 3,a_band{1});
I = speye(size(Lc{1}));
L = Lc{1};
M = E1*L - 4/dx^2*(I-E3);
% if dt=0.25*dx^2, then norm(A^k,inf) seems to grow all the time, might be
% unstable?
dt = 0.25*dx^2;
A = I + dt*M;


% kk = 1:1000;
% norm_Ak = zeros(size(kk));
% Ak = I;
% for i = kk
%     Ak1 = A*Ak;
%     tmp = (sum(abs(Ak1),2) - sum(abs(Ak),2));
%     %dAk = norm(Ak1-Ak,inf);
%     %dAk = max(max(Ak1)) - max(max(Ak));
%     dAk = min(min(Ak1)) - min(min(Ak));
%     Ak = Ak1;
%     norm_Ak(i) = norm(Ak,inf);
%     if mod(i,5) == 0
%         disp(['i=',num2str(i),'; norm(A^k,inf)=',num2str(norm_Ak(i)), ...
%               '; abs row sum diff: ', num2str(max(tmp)),'  ', num2str(min(tmp)), ...
%               '; diff: ', num2str(dAk)]);
%     end
% end