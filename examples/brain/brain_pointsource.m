%% Surface diffusion with reaction term and point sources
% on the brain
% u_t = lap_s u - alpha u /(1+u) + p.s.

% can speed up later runs if you set these to false
loaddata = true; build_matrices = true;
%loaddata = false; build_matrices = false;

if loaddata
dx = 0.05;

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


% parameters for brain reaction-diffusion
alpha = 0.1;   % coefficient of reaction term
pslam = 1;   % coefficient of point-sources forcing

% forcing function - gaussians around point sources
% first randomly choose some points, then find their closest points (this might bias eg. high curvature pts)
nsrcs = 3;
srcs = randi(length(band), nsrcs, 1);

% make a sum of gaussians
gsum = 0;
varsq = dx; % scale somehow
for srccount = 1:nsrcs
  si = srcs(srccount);
  gdist = (xg-cpxg(si)).^2 + (yg-cpyg(si)).^2 + (zg-cpzg(si)).^2;
  gsum = gsum + exp( -gdist/(2*varsq));
end
% cp-ext
gsum = E3*gsum;   % E1 or E3 here?

% viz this:
figure(2); clf;
set(0, 'CurrentFigure', 2);
clf;
gplot = Eplot*gsum;
trisurf(Faces,xp,yp,zp, gplot);
xlabel('x'); ylabel('y'); zlabel('z');
title( 'point sources' );
axis equal
shading interp
camlight left
lighting gouraud
colorbar
%fname = sprintf('frame_forcing.png');
%print('-dpng', fname)
drawnow()



%% Do calculation
%u0 = zeros(size(xg));
u0 = ones(size(xg));
u = u0;



Tf = 15000;
dt = 0.2*dx^2;
%dt = 1;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps


for kt = 1:numtimesteps
  % explicit Euler timestepping

  %unew = u + dt*E1*(L*u) - dt*alpha*u./(1+u) - dt*pslam*(u-gsum) - dt*lambda*(u - E3*u);
  unew = u + dt*E1*(L*u) - dt*alpha*u./(1+u) + dt*pslam*gsum - dt*lambda*(u - E3*u);

  u = unew;
 
  % closest point extension
  %tic
  %u = E3*unew;
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
    %fname = sprintf('frame%.4d.png', kt);
    %print('-dpng', fname)
    drawnow()
  end
end
