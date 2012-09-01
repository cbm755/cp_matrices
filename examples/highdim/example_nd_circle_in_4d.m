
dx = 0.1;

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;
w1d = x1d;

nx = length(x1d);
ny = length(y1d);
nz = length(z1d);
nw = length(w1d);

[xx yy zz ww] = ndgrid(x1d, y1d, z1d, w1d);

% closest point to a circle in the xx,yy plane embedded in 4D.
[th, r] = cart2pol(xx, yy);
[cpx, cpy] = pol2cart(th, 1);
ZC = dx/3;
WC = dx/4;
cpz = ZC*ones(size(cpx));
cpw = WC*ones(size(cpx));
dist = sqrt((xx-cpx).^2 + (yy-cpy).^2 + (zz-cpz).^2 + (ww-cpw).^2);
u0 = cos(th);


% banding
dim = 4;
p = 3;       % max interpolation order
stenrad = 1; % max stencil radius for finite differences
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((stenrad+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpx = cpx(band); cpy = cpy(band); cpz = cpz(band); cpw = cpw(band);
xg = xx(band); yg = yy(band); zg = zz(band); wg = ww(band);
dist = dist(band);
u0 = u0(band);


E1 = interpn_matrix({x1d y1d z1d w1d}, [cpx cpy cpz cpw], 1, band);
E3 = interpn_matrix({x1d y1d z1d w1d}, [cpx cpy cpz cpw], p, band);

L = laplacian_nd_matrix({x1d y1d z1d w1d}, 2, band);


%% Construct an interpolation matrix for plotting on circle

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,100)';
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:);  yp = yp(:);
zp = ZC*ones(size(xp));
wp = WC*ones(size(xp));

Eplot = interpn_matrix({x1d, y1d, z1d, w1d}, [xp yp zp wp], p, band);

figure(1); clf;
figure(2); clf;
figure(3); clf;


%% Time-stepping for the heat equation
Tf = 2;
%dt = 0.2*dx^2;
dt = 0.5*dx;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

u = u0;


lambda = 2*dim/dx^2;
I = speye(size(L));
M1 = E1*L - lambda*(I - E3);
A = I - dt*M1;


for kt = 1:numtimesteps
  % explicit Euler timestepping
  %unew = u + dt*(L*u);
  %u = E3*unew;

  % MOL, penalty method
  unew = A \ u;
  u = unew;

  % implicit Euler with bilinear interp p=1
  %unew = u + dt*((E1*L)*unew);
  %unew = A \ u;

  % closest point extension with p=3
  %u = E1*unew;
  %u = unew;

  t = kt*dt;

  if ( (kt < 5) || (mod(kt,numtimesteps/20) == 0) || (kt == numtimesteps) )
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
    legend('numerical soln', 'exact soln', 'Location', 'SouthEast');
    error_circ_inf = max(abs( exactplot - circplot ))

    set(0, 'CurrentFigure', 3);
    clf;
    plot(thetas, circplot - exactplot);
    title( ['error at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('error');

    %pause(0);
    drawnow();
  end
end
