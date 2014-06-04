%% Biharmonic equation on a circle
% This example solves the biharmonic equation on a 2D circle, with
% initial condition $u = \cos(2\theta)$, and exact solution
% $u(t) = \cos(2\theta) \exp(-16t)$.
%
% Author: Jun Lu, 2013.


%% Construct a grid in the embedding space

% Define grid size dx. Smaller dx values can be chosen if we wish to
% refine the grid mesh to obtain better approximations.
dx=0.1;
% Make vectors of x, y, which denotes positions of x and y
% coordinates of the grid points.
x1d = (-2.0:dx:2.0)';
y1d = x1d;
nx = length(x1d);
ny = length(y1d);


%% Find closest points on the surface
% For each grid point (x,y), we store the closest point on the circle.
% whose position expressed as (cpxx,cpyy).

% Meshgrid is only needed for finding the closest points, not afterwards
[xx yy] = meshgrid(x1d, y1d);
% Function 'cpCircle' is used for finding the closest points on a circle.
[cpxx, cpyy, dist] = cpCircle(xx,yy);
% Use cpxxv and cpxxy to denote the vector of x,y coordinates of closest
% points of all grid points
cpxxv = cpxx(:); cpyyv = cpyy(:);


%% Banding: do calculation in a narrow band around the circle

% Specify the dimension
dim = 2;
% Determine the degree of interpolation used in calculating the bandwidth.
% p, q are degree of interpolating polynomials used in constructing the
% interpolation matrices Ep and Eq.
% Degree of interpolating polynomial in computing the bandwidth is chosen
% to be the maximum value of p and q.
p=5; q=5;
interpdeg = max(p,q);
% Determine the order of Laplacian, typically degree 2 or 4.
Lorder = 4;

% "band" is a vector of the indices of the points in the computation band.
% The formula for bw is found in [Ruuth & Merriman 2008] and the 1.0001 is
% a safety factor.

bw = 1.0001*sqrt((dim-1)*((interpdeg+1)/2)^2 + ((Lorder/2+(interpdeg+1)/2)^2));
band = find(abs(dist)<=bw*dx);

% Select all the closest points in the band and store their coordinate
% information in vectors 'cpxxb' and 'cpyyb'.
cpxxb = cpxxv(band); cpyyb = cpyyv(band);
% Store coordinate information of grid points inside the interpolation band
% in vectors 'xxb' and 'yyb'.
xxb = xx(band); yyb = yy(band);


%% Function u in the embedding space
% u is a function defined on the grid that solves the biharmonic problem.
% Map points on meshgrid onto plane polar coordinates.
% Specify the initial condition, here in terms of the
% parameterization $\theta$.
[th,r] = cart2pol(xxb,yyb);
u = cos(2*th);
initialu = u;


%% Construct an interpolation and differentiation matrices
% Interpolation matrices denoted by E interpolates data from the
% uniform grid 'x1d' x 'y1d', onto the points 'cpxxb', 'cpyyb', the
% coordinates of all closest points in the band.
Ep = interp2_matrix(x1d,y1d,cpxxb,cpyyb,p,band);
Eq = interp2_matrix(x1d,y1d,cpxxb,cpyyb,q,band);
L = laplacian_2d_matrix(x1d,y1d,Lorder,band);


%%
% Matrices Ep, Eq are used in the following two steps:
% forward euler time-stepping: unew = u_n - dt*L*Ep*L*u_n
% closest point extension: u = Eq*unew;


%% Construct an interpolation matrix for plotting on circle

% Plot grid on circle and use theta as the parameter:
thetas = linspace(0,2*pi,100)';
r = ones(size(thetas));
% Plot grid in Cartesian coordinates
[xp,yp] = pol2cart(thetas,r);
xpv = xp(:); ypv = yp(:);
Eplot = interp2_matrix(x1d,y1d,xpv,ypv,interpdeg,band);


%% Time-stepping for the biharmonic equation
% First, fetermine the stopping time and size of each time step
Tf=0.01;
dt=0.04*dx^4;
% Adjust for integer number of steps
tsteps=ceil(Tf/dt);
dt = Tf / tsteps;


%%
% Loop in time, occasionally plotting solutions and errors.
figure(1); clf;
figure(2); clf;
figure(3); clf;

for kt = 1:tsteps
  % explicit Euler time-stepping
  unew = u - dt*(L*(Ep*(L*u)));
  % closest point extension
  u = Eq*unew;
  t = kt*dt;

  % ploting code
  if ( (kt < 10) || (mod(kt,1000) == 0) || (kt == tsteps) )
    % plot over computation band
    set(0, 'CurrentFigure', 1); clf
    plot2d_compdomain(u,xxb,yyb,dx,dx,1);
    title( ['embedded domain: soln at time ' num2str(t) ...
            ', timestep #' num2str(kt)] );
    xlabel('x'); ylabel('y');
    %hold on
    plot(xpv,ypv,'k-', 'linewidth',2);
    axis equal;  axis tight;
    hold off

    % plot value on circle
    set(0,'CurrentFigure', 2); clf;
    circplot = Eplot*u;
    % NB: circplot is computed, which is an approximated solution on the
    % cicle rather than in the band, so that error can be computed by
    % comparing with the exact solution.
    exactplot = cos(2*thetas)*exp(-16*t);
    plot(thetas,circplot);
    title( ['soln at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('u');
    hold on;
    % plot analytic result
    plot(thetas, exactplot, 'r--');
    legend('explicit Euler', 'exact answer', 'Location', 'SouthEast');
    err = norm(exactplot-circplot,inf)

    set(0, 'CurrentFigure', 3); clf;
    plot(thetas,exactplot-circplot);
    title( ['error at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('error');
    drawnow();
  end
end

