%% Advection equation on a circle
% This example solves the advection equation on a 2D circle using a
% specified velocity field tangent to the circle.
%
%% Advection on a curve
% The surface PDE is
%
% $$ u_t + \vec{W} \cdot \nabla_S u = 0, $$
%
% or more generally
%
% $$ u_t + \textrm{div}_S (u \vec{W} ) = 0, $$
%
% where $\vec{W}$ is a tangent vector field.
%

%% Build a meshgrid
% Construct a grid in the embedding space

dx = 0.05/2;   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;

% meshgrid is only needed for finding the closest points, not afterwards
[x y] = meshgrid(x1d, y1d);


%% Find closest points on the surface
% For each point (x,y), we store the closest point on the circle
% (cpx,cpy)

[cpx, cpy, dist] = cpCircle(x, y);


%% Banding: do calculation in a narrow band around the circle
dim = 2;  % dimension
p = 3;    % interpolation order
bw = rm_bandwidth(dim, p);
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpx = cpx(band); cpy = cpy(band);
x = x(band); y = y(band);


%% Function u in the embedding space
% u is a function defined on the grid and because we know the the
% surface is a circle, we use the parameterization to define an
% initial condition.  We define a function with the exact solution
% for calculating error later.
[th, r] = cart2pol(cpx, cpy);
uex = @(th,t) (cos(th-t)).^3;
u0 = uex(th, 0);

%% Construct interpolation and differentiation matrices
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

E = interp2_matrix(x1d, y1d, cpx, cpy, p, band);

L = laplacian_2d_matrix(x1d, y1d, 2, band);
[Dxb,Dxf, Dyb,Dyf] = firstderiv_upw1_2d_matrices(x1d, y1d, band);
[Dxc,Dyc] = firstderiv_cen2_2d_matrices(x1d, y1d, band);

%% Velocity vector in embedded space
% Unlike the instrinsic surface diffusion equation $u_t = \Delta_S u$,
% the advection equation on a curve needs a velocity field: a least a
% speed in the tangent direction is required as part of the *model*.
% We then need that velocity field embedded and extended into the 2D
% domain.

w1 = -y;
w2 = x;
w1 = E*w1;
w2 = E*w2;

%% Construct an interpolation matrix for plotting
% plotting grid on circle, using theta as a parameterization

thetas = linspace(0,2*pi,512)';
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix(x1d, y1d, xp, yp, p, band);




%% Time-stepping

Tf = 5;
dt = 0.25*dx;
numsteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numsteps

figure(1); clf;  figure(2); clf;  figure(3); clf;

u = u0;
for kt = 1:numsteps
  %% more generally, could compute velocities here:
  % w1 = ...
  % w2 = ...
  % w1 = E*w1;
  % w2 = E*w2;

  %% upwinding discretization of u_t + div_S . (u*vec{w}) = 0
  rhs = -( (w1 < 0) .* (Dxf*(u.*w1)) + (w1 >= 0) .* (Dxb*(u.*w1)) + ...
           (w2 < 0) .* (Dyf*(u.*w2)) + (w2 >= 0) .* (Dyb*(u.*w2)) );

  %% or centered diff: likely unstable here without diffusion
  %rhs = -( Dxc*(u.*w1) + (Dyc*(u.*w2)) );

  % explicit Euler timestepping
  unew = u + dt*rhs;

  % closest point extension
  u = E*unew;

  t = kt*dt;

  if ( (kt < 10) || (mod(kt,10) == 0) || (kt == numsteps) )
    % plot over computation band
    plot2d_compdomain(u, x, y, dx, dx, 1)
    title(sprintf('embedded domain: soln at time %g, timestep #%d', t, kt));
    xlabel('x');  ylabel('y');
    plot(xp, yp, 'k-', 'linewidth', 2);

    % plot value on circle
    set(0, 'CurrentFigure', 2); clf;
    circplot = Eplot*u;
    exactplot = uex(thetas, t);
    plot(thetas, circplot);
    hold on;
    plot(thetas, exactplot, 'r--');
    title(['soln at time ' num2str(t) ' on circle']);
    xlabel('\theta'); ylabel('u');
    legend('explicit Euler', 'exact answer', 'Location', 'NorthEast');
    error_circ_inf = max(abs( exactplot - circplot ));

    set(0, 'CurrentFigure', 3); clf;
    plot(thetas, circplot - exactplot);
    title(['error at time ' num2str(t) ' on circle']);
    xlabel('\theta'); ylabel('error');

    fprintf('step %d of %d, max_err=%g\n', kt, numsteps, error_circ_inf);
    drawnow();
  end
end
