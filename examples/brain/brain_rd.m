%% Gray--Scott reaction-diffusion

% can speed up later runs if you set these to false
loaddata = true; build_matrices = true;
%loaddata = false; build_matrices = false;

if loaddata
dx = 0.025;

x1d=(-2.0:dx:2.0)';
y1d=x1d;
z1d=x1d;
nx=length(x1d);
ny=length(y1d);
nz=length(z1d);
% cp_matrices assumes something like the following for ordering, but
% its not necessary to actually build the meshgrid:
%[x3d,y3d,z3d]=meshgrid(x1d,y1d,z1d);
%[cpx,cpy,cpz,dist] = cpSphere(x3d,y3d,z3d);


%PlyFile = 'bunny.ply';
PlyFile = 'brain-lh_scale_0.ply';
%PlyFile = 'annies_pig.ply';
%PlyFile = 'bumpy_torus_scaled.ply';
disp( ['reading triangulation from "' PlyFile '"'] );
tic
[Faces, Vertices] = plyread(PlyFile, 'tri');
toc

disp('converting to closest point representation');
[IJK,DIST,CP,XYZ] = tri2cp(Faces, Vertices, dx, x1d(1));
i = IJK(:,1);
j = IJK(:,2);
k = IJK(:,3);
dist = DIST;
cpxg = CP(:,1);
cpyg = CP(:,2);
cpzg = CP(:,3);
xg = XYZ(:,1);
yg = XYZ(:,2);
zg = XYZ(:,3);



dim = 3;
p = 3;  % degree interp
order = 2;  % laplacian order, griddata hardcoded for 2.


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

end


if (build_matrices)
  %% discrete operators
  disp('building laplacian and interp matrices');
  E3 = interp3_matrix(x1d,y1d,z1d, cpxg, cpyg, cpzg, p, band);
  E1 = interp3_matrix(x1d,y1d,z1d, cpxg, cpyg, cpzg, 1, band);
  L = laplacian_3d_matrix(x1d,y1d,z1d, order, band,band);
  % iCPM matrix
  %M = lapsharp(L,E3);
  lambda = 2*dim / dx^2;
  I = speye(size(E3));
  M = E1*L - lambda*(I-E3);


  %% plotting grid
  xp = Vertices(:,1);
  yp = Vertices(:,2);
  zp = Vertices(:,3);
  Eplot = interp3_matrix(x1d,y1d,z1d, xp,yp,zp, 3, band);
  %Eplot1 = interp3_matrix(x1d,y1d,z1d, xp,yp,zp, 1, band);
  %Eplot0 = interp3_matrix(x1d,y1d,z1d, xp,yp,zp, 0, band);
end



% parameters and functions for Gray--Scott
% 120 works with 0.025
F = 0.054;  kk = 0.063;  nuu = 1/(3/dx)^2;  nuv = nuu/2;
%F = 0.054;  k = 0.063;  nuu = 1/120^2;  nuv = nuu/2;
f = @(u,v) (-u.*v.*v  +  F*(1-u));
g = @(u,v) ( u.*v.*v  -  (F+kk)*v);


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
%dt = 1;
dt = .2* (1/max(nuu,nuv)) * dx^2
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

figure(2); clf;

for kt = 1:numtimesteps
  % explicit Euler timestepping
  %tic
  unew = u + dt*f(u,v) + dt*nuu*(E1*(L*u)) - dt*nuu*lambda*(u - E3*u);
  vnew = v + dt*g(u,v) + dt*nuv*(E1*(L*v)) - dt*nuv*lambda*(v - E3*v);
  %toc

  u = unew;
  v = vnew;
  % closest point extension
  %tic
  %u = E3*unew;
  %v = E3*vnew;
  %toc

  t = kt*dt;

  if ( (mod(kt,50)==0) || (kt<=5) || (kt==numtimesteps) )
    disp([kt t]);
    set(0, 'CurrentFigure', 2);
    clf;
    uplot = Eplot*u;
    trisurf(Faces,xp,yp,zp, uplot);
    xlabel('x'); ylabel('y'); zlabel('z');
    title( ['u at time ' num2str(t)] );
    axis equal
    shading interp
    camlight left
    lighting gouraud
    colorbar
    fname = sprintf('frame%.4d.png', kt);
    print('-dpng', fname)
    drawnow()
  end
end
