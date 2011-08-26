%% Gray--Scott reaction-diffusion on the Stanford Bunny
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.


% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');
addpath('../surfaces/tri_bunny');
addpath('../surfaces/readply');

loaddata = 1;

if (loaddata == 1)
% NOTE: griddata needs to be generated for each choice of dx (by the
% tri2cp C code).  its loaded below.  dim, p, order, x1d, y1d, z1d
% must match this code too.
dx = 0.025;

x1d=(-2.0):dx:(2.0)';
y1d=x1d;
z1d=x1d;
nx=length(x1d);
ny=length(y1d);
nz=length(z1d);
% the griddata file assumes something like this for ordering, but
% its not necessary to actually do the meshgrid.
%[x3d,y3d,z3d]=meshgrid(x1d,y1d,z1d);
%[cpx,cpy,cpz,dist] = cpSphere(x3d,y3d,z3d);

dim = 3;
p = 3;  % degree interp
order = 2;  % laplacian order, griddata hardcoded for 2.

% ply file contains the triangles (for plotting, see below)
PlyFile = 'bunny_input.ply';
% griddata: processed from ply file by tri2cp (C code).  Contains
% only points in the narrow band.
GD = load( ['bunny_griddata_p' num2str(p) ...
           '_dx' num2str(dx) '.txt'] );
% plus one b/c they're C indices
i = GD(:,1) + 1;
j = GD(:,2) + 1;
k = GD(:,3) + 1;
dist = GD(:,4);
cpxg = GD(:,5);
cpyg = GD(:,6);
cpzg = GD(:,7);
xg = GD(:,8);
yg = GD(:,9);
zg = GD(:,10);


%%sanity checks
xtest = x1d(1) + (i-1)*dx;
ytest = y1d(1) + (j-1)*dx;
ztest = z1d(1) + (k-1)*dx;
if ( (norm(xtest - xg, inf) > 1e-14) |  ...
     (norm(ytest - yg, inf) > 1e-14) | ...
     (norm(ztest - zg, inf) > 1e-14) )
  error('sanity fail');
end
%[temp,N] = size(GD);
%for c=1:N

%end

% here is one place where meshgrid comes in: note ordering here.
band = sub2ind([ny,nx,nz], j,i,k);


build_matrices = true;
if (build_matrices)
  %% discrete operators
  disp('building laplacian and interp matrices');
  L = laplacian_3d_matrix(x1d,y1d,z1d, order, band,band);
  E = interp3_matrix_band(x1d,y1d,z1d, cpxg, cpyg, cpzg, p, band);
  % iCPM matrix
  M = lapsharp(L,E);

  %% plotting grid
  [Faces, Vertices] = plyread(PlyFile, 'tri');
  xp = Vertices(:,1);
  yp = Vertices(:,2);
  zp = Vertices(:,3);
  Eplot = interp3_matrix_band(x1d,y1d,z1d, xp,yp,zp, 3, band);
  Eplot1 = interp3_matrix_band(x1d,y1d,z1d, xp,yp,zp, 1, band);
  %Eplot0 = interp3_matrix_band(x1d,y1d,z1d, xp,yp,zp, 0, band);
end


end

% parameters and functions for Gray--Scott
% 120 works with 0.025
F = 0.054;  k = 0.063;  nuu = 1/(3/dx)^2;  nuv = nuu/2;
f = @(u,v) (-u.*v.*v  +  F*(1-u));
g = @(u,v) ( u.*v.*v  -  (F+k)*v);


%% Do some calculation
% Here its eigenvalues/eigenfunctions but could be timestepping
% or whatever using the cp_matrices.

%pert = (yg - min(yg)) / (max(yg) - min(yg));
pert = (xg - min(xg)) / (max(xg) - min(xg));
pert = pert + 1.0*(rand(size(xg))-0.5);
u0 = 1 - pert;
v0 = 0 + 0.5*pert;
u = u0;
v = v0;

Tf = 15000;
%dt = 0.2*dx^2;
dt = 1;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

for kt = 1:numtimesteps
    % explicit Euler timestepping
    tic
    unew = u + dt*f(u,v) + dt*nuu*(L*u);
    vnew = v + dt*g(u,v) + dt*nuv*(L*v);
    toc
    
    % closest point extension
    tic
    u = E*unew;
    v = E*vnew;
    toc
    
    t = kt*dt;
    
    if ( (mod(kt,25)==0) | (kt<=10) | (kt==numtimesteps) )
    % plot value on sphere
    set(0, 'CurrentFigure', 2);
    clf;
    uplot = Eplot*u;
    % swap y and z here
    trisurf(Faces,xp,-zp,yp, uplot);
    xlabel('x'); ylabel('y'); zlabel('z');
    title( ['u at time ' num2str(t)] );
    axis equal
    shading interp
    camlight left
    colorbar
    pause(0.01);
    end

end
