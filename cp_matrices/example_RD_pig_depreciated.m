%% Gray--Scott reaction-diffusion on a triangulated pig.
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

warning('use example_RD_tri.m instead');

% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');
addpath('../surfaces/tri_pig');
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
PlyFile = 'pig_loop2.ply';
% griddata: processed from ply file by tri2cp (C code).  Contains
% only points in the narrow band.
GD = load( ['~/svn/closestpoint/surfaces/tri_pig/pig_loop2_griddata_p' num2str(p) ...
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
if ( (norm(xtest - xg, inf) > 1e-14) ||  ...
     (norm(ytest - yg, inf) > 1e-14) || ...
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
  tic;
  L = laplacian_3d_matrix_test(x1d,y1d,z1d, order, band,band);
  toc;
  tic;
  E = interp3_matrix_test(x1d,y1d,z1d, cpxg, cpyg, cpzg, p);
  E = E(:,band);
  toc;
  % iCPM matrix
  M = lapsharp(L,E);

  %% plotting grid
  [Faces, Vertices] = plyread(PlyFile, 'tri');
  xp = Vertices(:,1);
  yp = Vertices(:,2);
  zp = Vertices(:,3);
  Eplot = interp3_matrix_test(x1d,y1d,z1d, xp,yp,zp, 3);
  Eplot = Eplot(:,band);
  Eplot1 = interp3_matrix_test(x1d,y1d,z1d, xp,yp,zp, 1);
  Eplot1 = Eplot1(:,band);
  %Eplot0 = interp3_matrix_band(x1d,y1d,z1d, xp,yp,zp, 0, band);
end


end

% parameters and functions for Gray--Scott
% 120 works with 0.025
F = 0.054;  k = 0.063;  nuu = 1/120^2;  nuv = nuu/2;
f = @(u,v) (-u.*v.*v  +  F*(1-u));
g = @(u,v) ( u.*v.*v  -  (F+k)*v);


%% Do some calculation
% Here its eigenvalues/eigenfunctions but could be timestepping
% or whatever using the cp_matrices.

[th,r,temp] = cart2pol(xg,yg,zg);
%pert = 0.5*exp(-(6*(zg-.1)).^2) + 0.5*rand(size(xg));
pert = 1*exp(-(6*(zg-0.05*cos(6*th))).^2); % + 0.5*rand(size(xg));
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
    
    if ( (mod(kt,25)==0) || (kt<=10) || (kt==numtimesteps) )
    % plot value on sphere
    figure(2); clf;
    uplot = Eplot*u;
    trisurf(Faces,xp,yp,zp, uplot);
    xlabel('x'); ylabel('y'); zlabel('z');
    title( ['u at time ' num2str(t)] );
    axis equal
    shading interp
    camlight left
    colorbar
    pause(0);
    end

end
