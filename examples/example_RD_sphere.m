%% Reaction-diffusion on a sphere

cpf = @cpSphere;
paramf = @paramSphere;
%cpf = @cpEllipsoid;
%paramf = @paramEllipsoid;

loaddata = 1;

if (loaddata == 1)
  dx = 0.05;      % grid size

  % make vectors of x, y, z positions of the grid
  x1d = (-2.0:dx:2.0)';
  y1d = x1d;
  z1d = x1d;
  nx = length(x1d);
  ny = length(y1d);
  nz = length(z1d);

  % meshgrid is only needed for finding the closest points, not afterwards
  [xx yy zz] = meshgrid(x1d, y1d, z1d);

  [cpx, cpy, cpz, dist] = cpf(xx,yy,zz);
  cpx = cpx(:); cpy = cpy(:); cpz = cpz(:);

  %% Banding: do calculation in a narrow band around the surface
  dim = 3;  % dimension
  p = 3;    % interpolation order
  % "band" is a vector of the indices of the points in the computation
  % band.  The formula for bw is found in [Ruuth & Merriman 2008] and
  % the 1.0001 is a safety factor.
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(abs(dist) <= bw*dx);

  % store closest points in the band;
  cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
  x = xx(band); y = yy(band); z = zz(band);


  %% discrete operators
  disp('building laplacian and interp matrices');
  L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band,band);
  E = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p, band);
  I = speye(size(E));

  %% plotting grid
  [xp,yp,zp] = paramf(256);
  % Eplot is a matrix which interpolations data onto the plotting grid
  Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), p, band);
end

figure(1); clf;

% u_t = f(u,g) + nuu*Lap u
% v_t = g(u,g) + nuv*Lap u

% parameters and functions for Gray--Scott
FF = 0.054;  kk = 0.063;  nuu = 1/(3/dx)^2;  nuv = nuu/3;
f = @(u,v) (-u.*v.*v  +  FF*(1-u));
g = @(u,v) ( u.*v.*v  -  (FF+kk)*v);


%% initial conditions - small perturbation from steady state
pert = 0.5*exp(-(10*(z-.1)).^2) + 0.5*rand(size(x));
u0 = 1 - pert;  v0 = 0.5*pert;
u = u0;  v = v0;

Tf = 10000;
dt = .2 * (1/max(nuu,nuv)) * dx^2
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps


figure(1);
sphplot = Eplot*u;
sphplot = reshape(sphplot, size(xp));
Hplot = surf(xp, yp, zp, sphplot);
title('initial u')
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
view(-10, 60)
axis off;
shading interp
camlight left
colorbar

%% Method-of-lines approach
% See [vonGlehn/Macdonald/Maerz 2013]
%lambda = 6*max(nuu,nuv)/(dx^2);
%Au = nuu*(E*L) - lambda*(I-E);
%Av = nuv*(E*L) - lambda*(I-E);


for kt = 1:numtimesteps
  %% MOL: explicit Euler timestepping
  %unew = u + dt*( E*f(u,v) + Au*u );
  %vnew = v + dt*( E*g(u,v) + Av*v );
  %u = unew;
  %v = vnew;
  %% MOL: without precomputing matrices
  %rhsu = nuu*(L*u) + f(u,v);
  %rhsv = nuv*(L*v) + g(u,v);
  %unew = u + dt*( E*rhsu - lambda*(u - E*u) );
  %vnew = v + dt*( E*rhsv - lambda*(v - E*v) );
  %u = unew;
  %v = vnew;

  %% Ruuth-Merriman
  rhsu = nuu*(L*u) + f(u,v);
  rhsv = nuv*(L*v) + g(u,v);
  unew = u + dt*rhsu;
  vnew = v + dt*rhsv;
  u = E*unew;
  v = E*vnew;

  t = kt*dt;

  if ( (mod(kt,25)==0) || (kt<=10) || (kt==numtimesteps) )
    disp([kt t]);
    sphplot = Eplot*u;
    sphplot = reshape(sphplot, size(xp));
    set(0, 'CurrentFigure', 1);
    set(Hplot, 'CData', sphplot);
    title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    drawnow;
  end
end

