%% Eigenvalue problem on a triangulated pig.
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

% This example computes eigenvalues and eigenfunctions of the
% Laplace--Beltrami operator on a pig.  Demonstates how to load cp
% banded grid data from an external file


dx = 0.05;

% ply file contains the triangles
disp('reading plyread');
%PlyFile = 'bunny.ply';
PlyFile = 'pig_loop2.ply';
%PlyFile = 'annies_pig.ply';
[Faces, Vertices] = plyread(PlyFile, 'tri');

disp('running tri2cp');
[IJK,DIST,CP,XYZ] = tri2cp(Faces, Vertices, dx, -2);
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

x1d=-2.0:dx:2.0;
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
  L = laplacian_3d_matrix(x1d,y1d,z1d, order, band,band);
  E = interp3_matrix(x1d,y1d,z1d, cpxg, cpyg, cpzg, p, band);
  % iCPM matrix
  M = lapsharp(L,E);

  %% plotting grid
  xp = Vertices(:,1);
  yp = Vertices(:,2);
  zp = Vertices(:,3);
  Eplot = interp3_matrix(x1d,y1d,z1d, xp,yp,zp, 3, band);
  %Eplot1 = interp3_matrix(x1d,y1d,z1d, xp,yp,zp, 1, band);
  %Eplot0 = interp3_matrix(x1d,y1d,z1d, xp,yp,zp, 0, band);
end




%% Do some calculation
% Here its eigenvalues/eigenfunctions but could be timestepping
% or whatever using the cp_matrices.

tic
[V,D] = eigs(-M, 10, 'sm');
%[V,D] = eigs(-M, 6, 0.75);
evtime = toc
D = diag(D);
[Lambda,I] = sort(abs(D));
Lambda

figure(2); clf;
for i=1:4
  lambda = Lambda(i);
  eigenvec = V(:,I(i));

  uplot = Eplot*E*real(eigenvec);
  subplot(2,2,i);
  trisurf(Faces,xp,yp,zp, uplot);
  xlabel('x'); ylabel('y'); zlabel('z');
  title(['eigenvalue = ' num2str(lambda)]);
  axis equal
  shading interp
  camlight left
end

