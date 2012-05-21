%% Heat equation on a Yin Yang
% testing mixed codimension

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');

% add functions for finding the closest points
addpath('../surfaces');


%% Construct a grid in the embedding space

dx = 0.02;                   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;

nx = length(x1d);
ny = length(y1d);



%% Find closest points on the surface
% For each point (x,y), we store the closest point on the circle
% (cpx,cpy)

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy] = meshgrid(x1d, y1d);
[cpx, cpy, dist] = cpYinYang(xx,yy);
% make into vectors
cpxg = cpx(:); cpyg = cpy(:);


%% Banding: do calculation in a narrow band around the circle
dim = 2;  % dimension
p = 4;    % interpolation degree
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band);
xg = xx(band); yg = yy(band);


%% Function u in the embedding space
% u is a function defined on the grid (eg heat if solving the heat
% equation)
[thg, rg] = cart2pol(xg,yg);
u = cos(thg);
%u = (yg < 0) & (rg > 0.97)
initialu = u;       % store initial value


%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

E = interp2_matrix(x1d, y1d, cpxg, cpyg, p, band);

% e.g., closest point extension:
%u = E*u;

%% Create Laplacian matrix for heat equation

order = 2;  % Laplacian order: bw will need to increase if changed
L = laplacian_2d_matrix(x1d,y1d, order, band);


%% Construct an interpolation matrix for plotting on circle

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,500)';
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix(x1d, y1d, xp, yp, p, band);


figure(1); clf;
figure(2); clf;
%figure(3); clf;


%% Time-stepping for the heat equation
Tf = 2;
dt = 0.2*dx^2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

I = speye(size(L));
gamma = 4/dx^2;

for kt = 1:numtimesteps
  % explicit Euler timestepping
  unew = u + dt*(L*u);
  % closest point extension
  u = E*unew;

  %unew = u + dt*(E*(L*u)) - dt*gamma*(u - E*u);
  %u = unew;

  srcloc = -0.5*pi;
  II = find( (abs(thg - srcloc) < 0.02) & (rg > 0.9) );
  %II = find((xg == 0) & (cpyg == -1));
  u(II) = 1;
  t = kt*dt;

  if ( (kt < 10) || (mod(kt,500) == 0) || (kt == numtimesteps) )
    % plot over computation band
    plot2d_compdomain(u, xg, yg, dx, dx, 1)
    title( ['embedded domain: soln at time ' num2str(t) ...
            ', timestep #' num2str(kt)] );
    xlabel('x'); ylabel('y');
    %hold on
    plot(xp,yp,'k-', 'linewidth', 2);
    %axis equal;  axis tight;

    % plot value on circle
    set(0, 'CurrentFigure', 2);
    clf;
    circplot = Eplot*u;
    icplot = cos(thetas);
    plot(thetas, circplot);
    title( ['soln on outer circle at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('u');
    hold on;
    % plot analytic result
    plot(thetas, icplot, 'r--');
    legend('explicit Euler', 'initial condition', 'Location', 'SouthEast');
    %error_circ_inf = max(abs( exactplot - circplot ))
    plot((2*pi+srcloc)*[1 1], ylim, 'k--', 'linewidth', 2);
    plot((pi)*[1 1], ylim, 'k--');

    %set(0, 'CurrentFigure', 3);
    %clf;
    %plot(thetas, circplot - exactplot);
    %title( ['error at time ' num2str(t) ', on circle'] );
    %xlabel('theta'); ylabel('error');

    %pause(0);
    drawnow();
  end
end
