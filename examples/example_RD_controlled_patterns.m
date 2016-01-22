%% Gray--Scott reaction-diffusion on an ellipsoid
%  with a forcing term to control the patterns


% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');

% set up movie
% aviobj = avifile('RD_ellipsoid.avi', 'compression', 'None');

loaddata = 1;
build_matrices = 1;

if (loaddata == 1)

%%
% Construct a grid in the embedding space

dx = 0.025;      % grid size

% make vectors of x, y, z positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

nx = length(x1d);
ny = length(y1d);
nz = length(z1d);


%% Find closest points on the surface
% For each point (x,y,z), we store the closest point on the surface
% (cpx,cpy,cpz)

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy zz] = meshgrid(x1d, y1d, z1d);
% Using an ellipsoid
% function cpEllipsoid for finding the closest points
[cpx, cpy, cpz, dist] = cpEllipsoid(xx,yy,zz);
% make into vectors
cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);


%% Banding: do calculation in a narrow band around the surface
dim = 3;  % dimension
p = 3;    % interpolation order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
xg = xx(band); yg = yy(band); zg = zz(band);


% read in surface (if using a triangulated surface)
%PlyFile = 'bunny.ply';
%PlyFile = 'bumpy_torus_scaled.ply';
%disp( ['reading triangulation from "' PlyFile '"'] );
%[Faces, Vertices] = plyread(PlyFile, 'tri');


if (build_matrices)
  %% discrete operators
  disp('building laplacian and interp matrices');
  L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band,band);
  E = interp3_matrix(x1d,y1d,z1d, cpxg, cpyg, cpzg, p, band);
  I = speye(size(E));

  %% plotting grid
  % plotting grid on ellipsoid, based on parameterization
    [xp,yp,zp] = paramEllipsoid(64);
    xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
    % Eplot is a matrix which interpolations data onto the plotting grid
    Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);

    % plotting matrix at higher res
    [xph,yph,zph] = paramEllipsoid(128);
    xph1 = xph(:); yph1 = yph(:); zph1 = zph(:);
    EplotHigh = interp3_matrix(x1d, y1d, z1d, xph1, yph1, zph1, p, band);
    
end

end

figure(7); clf;

%% Grey-Scott with a forcing term

% u_t = f(u,g) + nuu*Lap u + lambda * Chi(x) * (u - 1);
% v_t = g(u,g) + nuv*Lap u + lambda * Chi(x) * (v - 0);

% forcing strength
gslam = -0.05/5;
%gslam = 0;

%% Create forcing function from loaded image
Im = imread('RD_mask.png');

%dmask = A(:,:,1);  % red channel for color;
dmask = Im;  % bw: use whole thing

% image has x,y swapped
[vertres,horzres] = size(dmask);

dmaskd = double(dmask);

% plane of image projection
surf_plane = 0;
  
%xmax = max(xg)
%xmin = min(xg)
%ymax = max(yg)
%ymin = min(yg)

% scale to surface size (change these depending on surface)
xmin =         -1.6;
xmax =          1.6;
ymin =         -1.2;
ymax =          1.2;

xi = round((horzres)*(xg - xmin) / (xmax - xmin));
yi = round((vertres)*(yg - ymin) / (ymax - ymin));
yi = vertres - yi;

chi = zeros(length(xg),1);


% create forcing function
for it = 1:length(xg)
    
  if (zg(it) >= surf_plane)     
    % front of the surface - image should appear here
    chi(it) = (255-dmaskd(yi(it),xi(it))) / 255.0;
    if (chi(it) > 0.9)
      chi(it) = 1; % words
    else
      chi(it) = 0; % not words
    end
    
  else
      % no forcing
      chi(it) = 0;
  end

end

% make sure forcing is a cp extension
chi = E*chi;
% re threshold
chi = chi>0.9;

% make second mask for changed word
%Im2 = imread('RD_mask2.png');

% if no change, use same image
Im2 = Im;

dmask = Im2;  % bw: use whole thing
% image has x,y swapped
[vertres,horzres] = size(dmask);

dmaskd = double(dmask);

xi = round((horzres)*(xg - xmin) / (xmax - xmin));
yi = round((vertres)*(yg - ymin) / (ymax - ymin));
yi = vertres - yi;

chi2 = zeros(length(xg),1);

% create second forcing function
for it = 1:length(xg)
    
  if (zg(it) >= surf_plane)     
    % front of the surface - image should appear here
    chi2(it) = (255-dmaskd(yi(it),xi(it))) / 255.0;
    if (chi2(it) > 0.9)
      chi2(it) = 1; % words
    else
      chi2(it) = 0; % not words
    end
    
  else
      % no forcing
      chi2(it) = 0;
  end

end

% make sure forcing function is a cp extension
chi2 = E*chi2;
% re threshold
chi2 = chi2>0.9;


%% plot forcing function
figure(10); clf;
forceplot = EplotHigh*(1-chi);
forceplot = reshape(forceplot, size(xph));
surf(xph, yph, zph, forceplot);
title('forcing function');
xlabel('x'); ylabel('y'); zlabel('z');
%caxis([-1.05 1.05]);   % lock color axis
axis equal; shading interp
camlight left; colorbar
% colormap(fireice) - a nice colormap
drawnow()


%%

% parameters and functions for Gray--Scott
% 120 works with 0.025
F = 0.054;  kk = 0.063;  nuu = 1/(3/dx)^2;  nuv = nuu/2; % parameters for
%sphere
%F = 0.054;  kk = 0.063;  nuu = 1/(2/dx)^2;  nuv = nuu/2;
%F = 0.054;  k = 0.063;  nuu = 1/120^2;  nuv = nuu/2;
f = @(u,v) (-u.*v.*v  +  F*(1-u));
g = @(u,v) ( u.*v.*v  -  (F+kk)*v);


%% initial conditions - small perturbation from steady state
[th,r,temp] = cart2pol(xg,yg,zg);
pert = 0.5*exp(-(10*(zg-.1)).^2) + 0.5*rand(size(xg));
%pert = 1*exp(-(6*(zg-0.05*cos(6*th))).^2);% + 0.5*rand(size(xg));
u0 = 1 - pert;
v0 = 0 + 0.5*pert;
u = u0;
v = v0;

Tf = 5000;
%dt = 0.2*dx^2;
%dt = 1;
dt = .1* (1/max(nuu,nuv)) * dx^2
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

%% cpmol

lambda = 4*nuu/(dx^2); % similar to (scaled) lapsharp
%lambda = 1/dt;

Au = nuu*(E*L) - lambda*(I-E);
Av = nuv*(E*L) - lambda*(I-E);

%%

figure(2); clf;

for kt = 1:numtimesteps
  % explicit Euler timestepping
  
  if (kt < 10) % no forcing
    unew = u + dt*( E*f(u,v) + Au*u);
    vnew = v + dt*( E*g(u,v) + Av*v);
  elseif (kt < 3000) % turn on first forcing function
    unew = u + dt*( E*f(u,v) + Au*u + gslam * chi.*(u-0.3));
    vnew = v + dt*( E*g(u,v) + Av*v + gslam * chi.*(v-0.6));
  else    % turn on second forcing function
    unew = u + dt*( E*f(u,v) + Au*u + gslam * chi2.*(u-0.3));
    vnew = v + dt*( E*g(u,v) + Av*v + gslam * chi2.*(v-0.6));
  end

  
 u = unew;
 v = vnew;

  t = kt*dt;

  if ( (mod(kt,20)==0) || (kt<=50) || (kt==numtimesteps) )
    disp([kt t]);
    
    % plot surface: u function
    set(0, 'CurrentFigure', 2);
    clf;
    sphplot = EplotHigh*u;
    sphplot = reshape(sphplot, size(xph));
    surf(xph, yph, zph, sphplot);
    title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    xlabel('x'); ylabel('y'); zlabel('z');
    %caxis([-1.05 1.05]);   % lock color axis
    axis equal
    view(-10, 60)
    axis off;
    shading interp
    camlight left
    %colormap(fireice);
    colorbar
    drawnow()
    
    
    % plot surface: v function
    set(0, 'CurrentFigure', 7);
    clf;
    sphplot = EplotHigh*v;
    sphplot = reshape(sphplot, size(xph));
    surf(xph, yph, zph, sphplot);
    %title( ['v at time ' num2str(t) ', kt= ' num2str(kt)] );
    xlabel('x'); ylabel('y'); zlabel('z');
    %caxis([-1.05 1.05]);   % lock color axis
    axis equal
    view(-10, 60)
    axis off;
    shading interp
    camlight left; colorbar
    %colormap(fireice);
    drawnow()
    
    % movie
    % F = getframe;
    % aviobj = addframe(aviobj,F);
  end
end


% close movie file
%aviobj = close(aviobj);
