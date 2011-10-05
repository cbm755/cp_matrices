%% Heat equation on a circle
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

% This example demonstrates two bands as in the implicit CP paper
% [Macdonald, Ruuth 2009]


%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');

% add functions for finding the closest points
addpath('../surfaces');


global ICPM2009BANDINGCHECKS

% this is a bit dangerous: will break other less tightly banded
% codes, turn it off later
ICPM2009BANDINGCHECKS = 1;

%%
% 2D example on a circle
% Construct a grid in the embedding space

dx = 0.25/2;   % grid size

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
% function cpCircle for finding the closest points on a circle
[cpx, cpy, dist] = cpCircle(xx,yy);
% make into vectors
%cpxg = cpx(:); cpyg = cpy(:);


%% Banding: do calculation in narrow bands
% We start by defining an "initial band" which we'll refine later.
% Typically, in practice, we have some scattered input which we'll
% assume satisfies the bandwidth formula below (e.g., output from
% tri2cp)
dim = 2;  % dimension
p = 3;    % interpolation order
order = 2;  % Laplacian order: bw will need to increase if changed
fd_stenrad = order/2;  % Finite difference stencil radius
% The formula for bw is found in [Ruuth & Merriman 2008] and the
% 1.0002 is a safety factor.
bw = 1.0002*sqrt((dim-1)*((p+1)/2)^2 + ((fd_stenrad+(p+1)/2)^2));
% start with a rough band and refine later, here just find the
% indicies of all points within bandwidth of the surface.
band_init = find(abs(dist) <= bw*dx);
% the corresponding closest points
cpxg_init = cpx(band_init); cpyg_init = cpy(band_init);
xg_init = xx(band_init); yg_init = yy(band_init);


%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.
disp('Constructing interpolation matrix');
% various alternatives, here we use the "full" E matrix...
Etemp = interp2_matrix(x1d, y1d, cpxg_init, cpyg_init, p);
tic; [i,j,S] = find(Etemp); toc;
tic; innerband = unique(j); toc;


%% Create Laplacian matrix for heat equation
% in general, want the biggest stencil here so the others fit too
Ltemp = laplacian_2d_matrix(x1d,y1d, order, innerband, band_init);


%% The outerband
% We find a narrow outerband by using the column-space of L.
tic; [i,j,S] = find(Ltemp); toc;
tic; outerbandtemp = unique(j); toc;
% indices are into band_init (the original columns of L),
% look them up to get the outerband in terms of the
% meshgrid(x1d,y1d) indices.
outerband = band_init(outerbandtemp);

cpxgout = cpxg_init(outerbandtemp); cpygout = cpyg_init(outerbandtemp);
xgout = xg_init(outerbandtemp); ygout = yg_init(outerbandtemp);

L = Ltemp(:, outerbandtemp);
E = Etemp(outerbandtemp, innerband);
clear Ltemp Etemp outerbandtemp  % optional, erase the originals


%% Could instead regenerate everything, now that we have the two bands
%cpxgin3 = cpx(innerband); cpygin3 = cpy(innerband);
%cpxgout3 = cpx(outerband); cpygout3 = cpy(outerband);
%E3 = interp2_matrix_band(x1d, y1d, cpxgout3, cpygout3, p, innerband);
%L3 = laplacian_2d_matrix(x1d,y1d, order, innerband, outerband);


%% Other operators
% Note: stencils need to be a subset of the one used to generate
% the outerband above
[Dxb,Dxf,Dyb,Dyf] = firstderiv_upw1_2d_matrices(x1d,y1d, innerband, outerband);
[Dxc,Dyc] = firstderiv_cen2_2d_matrices(x1d,y1d, innerband, outerband);
[Dxx,Dyy] = secondderiv_cen2_2d_matrices(x1d,y1d, innerband, outerband);
% e.g., this one needs diagonals and might not work
%Dxy = secondderiv_mixcen2_2d_matrix(x1d,y1d, innerband, outerband);

%% Restriction operator
% used to extract inner values from an outer band vector.  There
% is probably a slick loop-free way to do this.
innerInOuter = zeros(size(innerband));
R = sparse([],[],[],length(innerband),length(outerband),length(innerband));
for i=1:length(innerband)
  I = find(outerband == innerband(i));
  innerInOuter(i) = I;
  R(i,I) = 1;
end
% closest points of the inner band
cpxgin = R*cpxgout;  cpygin = R*cpygout;
xgin = R*xgout;  ygin = R*ygout;


%% Diagonal splitting for iCPM
M = lapsharp_unordered(L, E, R);



%% Construct an interpolation matrix for plotting on circle
% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,1000);
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix_band(x1d, y1d, xp, yp, p, innerband);


% after building matrices, don't need this set
ICPM2009BANDINGCHECKS = 0;


%% Function u in the embedding space, initial conditions
% u is a function defined on the grid (e.g. heat)
[thg, rg] = cart2pol(cpxgin,cpygin);
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
    plot2d_compdomain(u, xgin, ygin, dx, dy, 1);
    hold on;
    plot(xp,yp,'k-', 'linewidth', 2);
    title( ['embedded domain: soln at time ' num2str(t) ...
            ', timestep #' num2str(kt)] );

    % plot value on circle
    figure(2); clf;
    circplot = Eplot*u;
    plot(thetas, circplot);
    title( ['soln at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('u');
    hold on;
    % plot analytic result
    plot(thetas, exp(-t)*cos(thetas), 'r--');
    plot(thetas, Eplot*u0, 'g-.');
    legend('explicit Euler', 'exact answer', 'initial condition ', ...
           'Location', 'SouthEast');
    %error_circ_inf = max(abs( exp(-t)*cos(thetas) - circplot' ));
    error_circ_inf = max(abs( uexactfn(t,thetas) - circplot' ));
    [dx dt t error_circ_inf]

    pause(0);
  end
end
