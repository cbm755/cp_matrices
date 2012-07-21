%% Heat equation on a circle
% This example solves the heat equation on a 2D circle, with initial
% conditions u = cos(theta), and exact solution u(t) =
% exp(-t)*cos(theta)


%% Using cp_matrices

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');

% add functions for finding the closest points
addpath('../surfaces');


%% Construct a grid in the embedding space

dx = 0.1/2^2;                   % grid size

%0 0.022321
%1 0.0053188
%2 0.27585

%1 0.011468
%2 incons

%1 0.016548
%2 0.007989
%3 0.0066063
%4 0.0071652

% 2-proj
%2 0.006579
%3 0.0046597
%4 0.0052347

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
% function cpCircle for finding the closest points on a circle
[cpxx, cpyy, dist] = cpCircle(xx,yy);
% make into vectors
cpxx = cpxx(:); cpyy = cpyy(:);


%% Banding: do calculation in a narrow band around the circle
dim = 2;  % dimension
p = 3;    % interpolation degree
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) <= 2.5*bw*dx);  % just some overestimate
bw2 = 2*bw*dx;
band2 = find(abs(dist) <= bw2);


% store closest points in the band;
cpxg = cpxx(band); cpyg = cpyy(band);
xg = xx(band); yg = yy(band);


%% Function u in the embedding space
% u is a function defined on the grid (eg heat if solving the heat
% equation)
%[th, r] = cart2pol(xg,yg);
%u = cos(th);
%u0 = u;       % store initial value


%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

%E = interp2_matrix(x1d, y1d, cpxg, cpyg, p, band);
%E1 = interp2_matrix(x1d, y1d, cpxg, cpyg, 1, band);

% e.g., closest point extension:
%u = E*u;

%% Create Laplacian matrix for heat equation

order = 2;  % Laplacian order: bw will need to increase if changed

%L = laplacian_2d_matrix(x1d,y1d, order, band);


%% Construct an interpolation matrix for plotting on circle

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,100)';
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);
%Eplot = interp2_matrix(x1d, y1d, xp, yp, p, band);


figure(1); clf;
figure(2); clf;
figure(3); clf;
figure(4); clf;

%% Time-stepping for the heat equation
Tf = 1;
%dt = .25*dx^2;
dt = dx/2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps


%stop

%% split a function off for this
% given a big band (overestimate) and an inner band (possibly but not
% necessarily the minimal one), find the ghost points for a finite
% distance stencil
Ltemp = laplacian_2d_matrix(x1d,y1d, order, band2, band);
[i,j,aij] = find(Ltemp);
j2 = unique(j);
band2wghosts = band(j2);
band2ghosts = setdiff(band2wghosts, band2);

%% now reorder so ghost's are at the end
% maybe this function should return 2 or 3 cp_grids?
band3 = [band2; band2ghosts];
% TODO: not necessary to recalcalute, see ops_and_bands code
L = laplacian_2d_matrix(x1d,y1d, order, band2, band3);
assert(nnz(Ltemp)==nnz(L))
% TODO: xg, yg, cpxg, cpyg, etc are wrong.  cp_grid structure
% should be passed to this function and fixed here

% TODO: move later
%E = interp2_matrix(x1d, y1d, cpxx(band2), cpyy(band2), 3, band2);
%Eplot = interp2_matrix(x1d, y1d, xp, yp, 3, band2);
E3 = interp2_matrix(x1d, y1d, cpxx(band3), cpyy(band3), 3, band3);
E1 = interp2_matrix(x1d, y1d, cpxx(band3), cpyy(band3), 1, band3);
Eplot = interp2_matrix(x1d, y1d, xp, yp, 3, band3);


%% interpolation for the ghost points
% this would be the ghost cp_grid.  All these can be found in band2
% by assumption that band2 is at least the minimal grid
Eg1 = interp2_matrix(x1d, y1d, cpxx(band2ghosts), cpyy(band2ghosts), 3, band2);

% TODO: try an actual neumann condition
cpx_g = cpxx(band2ghosts);
cpy_g = cpyy(band2ghosts);
x_g = xx(band2ghosts);
y_g = yy(band2ghosts);
n1 = x_g - cpx_g;
n2 = y_g - cpy_g;
dist_g = dist(band2ghosts);

offset = abs(dist_g)-bw2;

assert(min(offset)>0)

plot2d_compdomain(x_g, x_g, y_g, dx, dx, 4);

if (1==1)
  x_p = x_g - 2*offset ./ abs(dist_g) .* n1;
  y_p = y_g - 2*offset ./ abs(dist_g) .* n2;
else
  nn1 = n1 ./ sqrt(n1.^2 + n2.^2);
  nn2 = n2 ./ sqrt(n1.^2 + n2.^2);
  x_p = x_g - sqrt(2)*dx .* nn1;
  y_p = y_g - sqrt(2)*dx .* nn2;
end
for i=1:length(band2ghosts)
  plot(x_p(i),y_p(i),'k.')
  plot([x_g(i) x_p(i)], [y_g(i) y_p(i)], 'k-')
end
Eg2 = interp2_matrix(x1d, y1d, x_p, y_p, 1, band3);
%Eg2 = interp2_matrix(x1d, y1d, cpx_g, cpy_g, 1, band3);
%pause
Ii = speye(size(L,1),size(L,2));
Ig = speye(length(band2ghosts));
Z = sparse(length(band2ghosts), length(band2));
B = [Z Ig];
B = B - Eg2;

A = [Ii - dt*L; B];

% store closest points in the band;
cpxg = cpxx(band3); cpyg = cpyy(band3);
xg = xx(band3); yg = yy(band3);

%% Function u in the embedding space
% u is a function defined on the grid (eg heat if solving the heat
% equation)
[th, r] = cart2pol(xg,yg);
u = cos(th);
u0 = u;       % store initial value



%I = speye(size(band2));
%A = I - dt*(E1*L);
%A = I - dt*(LB);


M = length(band2)+1;

%whband = u0;
%whband(M:end) = u0(M:end)+1;
%plot2d_compdomain(whband, xg, yg, dx, dx, 4)

for kt = 1:numtimesteps
  % explicit Euler timestepping
  if (1==0)
    v = u(1:M-1);
    vv = [v; Eg2*u];
    w = v + dt*(L*vv);

    Z = Eg2(:,M:end);
    nnz(Z);
    AA = Ig - Z;
    W = AA \ (Eg2(:,1:M-1)*w);
    unew = [w; W];
    u = unew;

  else % implicit Euler

    v = u;
    v(length(band2)+1:end) = 0;
    unew = A \ v;

    % closest point extension
    if (1==1)
      u = E1*unew;
    else
      u = unew;
    end
  end

  t = kt*dt;

  if ( (kt < 5) || (mod(kt,20) == 0) || (kt == numtimesteps) )
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
    exactplot = exp(-t)*cos(thetas);
    plot(thetas, circplot);
    title( ['soln at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('u');
    hold on;
    % plot analytic result
    plot(thetas, exactplot, 'r--');
    legend('explicit Euler', 'exact answer', 'Location', 'SouthEast');
    error_circ_inf = max(abs( exactplot - circplot ))

    set(0, 'CurrentFigure', 3);
    clf;
    plot(thetas, circplot - exactplot);
    title( ['error at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('error');

    %pause;
    drawnow();
  end
end
