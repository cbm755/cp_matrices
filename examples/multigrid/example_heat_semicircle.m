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


x0 = -3;
x1 = 3;
y0 = -3;
y1 = 3;

%%
% 2D example on a circle
% Construct a grid in the embedding space

%dx = 0.003125; % grid size
dx = 0.00625;
dx_coarsest = 0.2;   % coarsest grid size
x1d_coarsest = (x0:dx_coarsest:x1)';
y1d_coarsest = (y0:dx_coarsest:y1)';

dy = dx;
dt = dx/10;

dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed

bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

w = 1;
n1 = 2;
n2 = 1;

MAX = 1;

p_mg = 3;

cpf1 = @cpSemicircle;  paramf = @paramSemicircle;  cpf = @(x,y) cpbar_2d(x,y,cpf1);

%uexactfn = @(t,th) exp(-t)*sin(th) + exp(-9*t)*sin(3*th) + exp(-25*t)*sin(5*th);
uexactfn = @(t,th) exp(-t)*cos(th) + exp(-9*t)*cos(3*th);

% caution: rhsfn is of not use but for calling 'helper_set_variables' This
% is not good. Hope to modify.
rhsfn = @(th) uexactfn(0,th);

has_boundary = true;

%[a_x1d, a_y1d, a_xcp, a_ycp, a_band, Mc, Lc, Ec, V, F, A, a_bdyg] = ...
%       helper_set_variables(x0, x1, y0, y1, dx, dx_coarsest, dim, p, order, rhsfn, cpf, has_boundary);


[a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d] = ... 
    build_mg_cpgrid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary);
[Mc, Lc, Ec] = build_mg_cpmatrix(a_band, a_xcp, a_ycp, a_x1d, a_y1d, p, order);

F = [];
[Mc, Lc, Ec, F] = app_bnd(Mc, Lc, Ec, F, a_xcp, a_ycp, a_bdyg, 'dirichlet');

% For time-stepping Neumann boundary condition, DON'T do the following
% line of code
% [Mc, Lc, Ec, F] = app_bnd(Mc, Lc, Ec, F, a_xcp, a_ycp, a_bdyg, 'neumann');

n_level = length(a_band);
for i = 1:1:n_level
    Mc{i} = speye(size(Mc{i})) - dt*Mc{i};
    Lc{i} = speye(size(Lc{i})) - dt*Lc{i}; 
end 

[TMf2c, TMc2f] = helper_set_TM(a_x1d, a_y1d, a_xcp, a_ycp, a_band, a_bdyg, p, p_mg);

cpxg = a_xcp{1};
cpyg = a_ycp{1};
[thg, rg] = cart2pol(cpxg,cpyg);
u0 = uexactfn(0,thg);
u = u0;

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,pi,1000)';
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
Tf = 0.1;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;

gap = 200;
err_plot = zeros(numtimesteps,1);
err_plot1 = zeros(numtimesteps,1);
t_plot = dt:dt:Tf;

figure(1);
tic
for kt = 1:numtimesteps
    
    if kt == 1
        % the following line is to enforce dirichlet boundary condition at
        % the first time step.
        u = Ec{1}*u;
    end
    
    t = kt*dt;
    %unew = Mc{1} \ u;
    [unew flag] = gmres(Mc{1},u);
    
    u = unew;

    circplot = Eplot{1}*u;
    error_circ_inf = max(abs( uexactfn(t,thetas) - circplot )) / max(abs(uexactfn(t,thetas)));
    err_plot1(kt) = error_circ_inf;
    if ( (kt < 5) || (mod(kt,gap) == 0) || (kt == numtimesteps) )
        [t error_circ_inf]
    end
    %plotting
%     if ( (kt < 5) || (mod(kt,gap) == 0) || (kt == numtimesteps) )
%       t
%       % plot in the embedded domain: shows the computational band
%       plot2d_compdomain(u, XCP{1}, YCP{1}, dx, dy, 1);
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
figure(3);
tic
for kt = 1:numtimesteps
    
    if kt == 1
        % the following line is to enforce dirichlet boundary condition at
        % the first time step.
        u = Ec{1}*u;
    end
    t = kt*dt;
    %% implicit backward euler
    % set the initial value of multigrid as the solution of last time step
    [V, F] = helper_set_V_F(a_band, u);
    unew = gmg(Mc, Lc, Ec, V, F, TMf2c, TMc2f, a_band, a_bdyg, n1, n2, 1, w, uexactfn(t,thetas), Eplot, MAX);

    u = unew;

    circplot = Eplot{1}*u;
    error_circ_inf = max(abs( uexactfn(t,thetas) - circplot )) / max(abs(uexactfn(t,thetas)));
    err_plot(kt) = error_circ_inf;
    if ( (kt < 5) || (mod(kt,gap) == 0) || (kt == numtimesteps) )
        [t error_circ_inf]
    end
    %plotting
%     if ( (kt < 5) || (mod(kt,gap) == 0) || (kt == numtimesteps) )
%       t
%       % plot in the embedded domain: shows the computational band
%       plot2d_compdomain(u, XCP{1}, YCP{1}, dx, dy, 1);
%       hold on;
%       plot(xp,yp,'k-', 'linewidth', 2);
%       title( ['embedded domain: soln at time ' num2str(t) ...
%             ', timestep #' num2str(kt)] );
% 
%       % plot value on circle
%       figure(4); clf;
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

figure(5)
plot(t_plot, err_plot1, t_plot,err_plot,'r')
xlabel('time')
ylabel('relative error')
%title(['relative error of e^{-t}cos(\theta)+e^{-9t}cos(3\theta), dx = ', num2str(dx)])
title(['relative error of e^{-t}sin(\theta)+e^{-9t}sin(3\theta), dx = ', num2str(dx)])
legend('matlab','multigrid')
