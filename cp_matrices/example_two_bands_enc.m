%% Heat equation on a circle
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

% This example demonstrates two bands as in the implicit CP paper,
% [Macdonald, Ruuth 2009] using "ops_and_bands2d" helper.


%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');

% add functions for finding the closest points
addpath('../surfaces');


%% 2D example on a circle
% Construct a grid in the embedding space

dx = 0.05;   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2:dx:2)';
y1d = (-1.4:dx:1.4)';
dy = dx;

nx = length(x1d);
ny = length(y1d);


%% Find closest points on the surface
% For each point (x,y), we store the closest point on the circle
% (cpx,cpy)

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy] = meshgrid(x1d, y1d);
% find the closest points
cpf = @cpCircle;  param = @paramCircle;
%cpf = @cpSemicircle;  param = @paramSemicircle;
[cpx2d, cpy2d, dist2d] = cpf(xx, yy);
%[cpx2d, cpy2d, dist2d, bdy2d] = cpbar_2d(xx, yy, cpf);


%% Banding: do calculation in narrow bands
% We start by defining an "initial band" which we'll refine later.
% Typically, in practice, we have some scattered input which we'll
% assume satisfies the bandwidth formula below (e.g., output from
% tri2cp)
dim = 2;  % dimension
p = 3;    % interpolation degree
order = 2;  % Laplacian order: bw will need to increase if changed
fd_stenrad = order/2;  % Finite difference stencil radius
% The formula for bw is found in [Ruuth & Merriman 2008] and the
% 1.0002 is a safety factor.
bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((fd_stenrad+(p+1)/2)^2));
% start with a rough band and refine later, here just find the
% indicies of all points within bandwidth of the surface.
bandinit = find(abs(dist2d) <= bw*dx);
% the corresponding closest points
cpxinit = cpx2d(bandinit); cpyinit = cpy2d(bandinit);
xinit = xx(bandinit); yinit = yy(bandinit);


%% Operators and banding
% take the initial band (which might be too wide) and generate the
% inner and outer band.  Also build Laplacian, Extension matrix and
% Restriction operator.

% this is a bit dangerous: will break other less tightly banded
% codes, turn it off later
global ICPM2009BANDINGCHECKS
ICPM2009BANDINGCHECKS = 1;

[L, E, R, iband, oband, iband2, oband2] = ...
    ops_and_bands2d(x1d,y1d, xinit,yinit, cpxinit,cpyinit, bandinit, p, order);

x = xinit(iband);
y = yinit(iband);
cpx = cpxinit(iband);
cpy = cpyinit(iband);
xout = xinit(oband);
yout = yinit(oband);
cpxout = cpxinit(oband);
cpyout = cpyinit(oband);
%bdy = bdyinit(iband);
%bdyout = bdyinit(oband);


%% Other operators
% If we need other operators and they will fit in the band, build them
% now.  (Note: stencils need to be a subset of the one used to
% generate the outerband above)
[Dxb,Dxf,Dyb,Dyf] = firstderiv_upw1_2d_matrices(x1d,y1d, iband2, oband2);
[Dxc,Dyc] = firstderiv_cen2_2d_matrices(x1d,y1d, iband2, oband2);
[Dxx,Dyy] = secondderiv_cen2_2d_matrices(x1d,y1d, iband2, oband2);
% e.g., this one needs diagonals and might not work (i.e., won't!)
%Dxy = secondderiv_mixcen2_2d_matrix(x1d,y1d, iband2, oband2);


%% Diagonal splitting for iCPM
M = lapsharp_unordered(L, E, R);


%% Construct an interpolation matrix for plotting on circle
% plotting grid on circle, using theta as a parameterization
%thetap = linspace(0,2*pi,1000)';
%r = ones(size(thetap));
% plotting grid in Cartesian coords
%[xp,yp] = pol2cart(thetap,r);
%xp = xp(:); yp = yp(:);
[xp,yp,thp] = param(600);
%Eplot = interp2_matrix_band(x1d, y1d, xp, yp, p, iband2);
Eplot = interp2_matrix(x1d, y1d, xp, yp, p, iband2);

% after building matrices, don't need this set
ICPM2009BANDINGCHECKS = 0;



%% Function u in the embedding space, initial conditions
% u is a function defined on the grid (e.g. heat)
[thg, rg] = cart2pol(cpx,cpy);
u0 = cos(thg);
uexactfn = @(t,th) exp(-t)*cos(th);
u = u0;

% for making full dense plots
%ufull = zeros(ny,nx);
%[th, r] = cart2pol(xx,yy);
%ufull = cos(th);


%% Time-stepping for the heat equation

implicit = 1
Tf = 0.25;
if implicit
  dt = dx/10;
else
  dt = 0.2*dx^2;
end
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps
if implicit
  % make time-stepping matrix
  I = speye(size(M));
  A = I - dt*M;
end

figure(1);

for kt = 1:numtimesteps
  if (implicit)
    %% implicit backward euler
    %unew - uold = dt*M*unew
    unew = A \ u;
  else
    %% explicit forward euler
    % closest point extension
    uext = E*u;
    % timestep
    unew = R*uext + dt*(L*uext);
    %unew = R*uext + dt*M*(R*uext);
  end

  u = unew;

  t = kt*dt;

  % plotting
  if ( (kt < 5) || (mod(kt,200) == 0) || (kt == numtimesteps) )
    % plot in the embedded domain: shows the computational band
    plot2d_compdomain(u, x, y, dx, dy, 1);
    hold on;
    plot(xp, yp, 'k-', 'linewidth', 2);
    title( ['embedded domain: soln at time ' num2str(t) ...
            ', timestep #' num2str(kt)] );

    % plot value on circle
    figure(2); clf;
    %subplot(2,1,2);
    circplot = Eplot*u;
    plot(thp, circplot);
    title( ['soln at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('u');
    hold on;
    % plot analytic result
    plot(thp, exp(-t)*cos(thp), 'r--');
    plot(thp, Eplot*u0, 'g-.');
    legend('explicit Euler', 'exact answer', 'initial condition ', ...
           'Location', 'SouthEast');
    %error_circ_inf = max(abs( exp(-t)*cos(thp) - circplot ));
    error_circ_inf = max(abs( uexactfn(t,thp) - circplot ));
    [dx dt t error_circ_inf]

    pause(0);
  end
end
