% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');

% NOTE: griddata needs to be generated for each choice of dx (by the
% tri2cp C code).  its loaded below.  dim, p, order, x1d, y1d, z1d
% must match this code too.
dx = 0.1;

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
GD = load( ['pig_loop2_griddata_p' num2str(p) ...
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
  M = diagSplit(L,E);

  %% plotting grid
  [Faces, Vertices] = plyread(PlyFile, 'tri');
  xp = Vertices(:,1);
  yp = Vertices(:,2);
  zp = Vertices(:,3);
  Eplot = interp3_matrix_band(x1d,y1d,z1d, xp,yp,zp, 3, band);
  Eplot1 = interp3_matrix_band(x1d,y1d,z1d, xp,yp,zp, 1, band);
  %Eplot0 = interp3_matrix_band(x1d,y1d,z1d, xp,yp,zp, 0, band);
end




%% Do some calculation
% Here its eigenvalues/eigenfunctions but could be timestepping
% or whatever using the cp_matrices.

tic
%[V,D] = eigs(-M, 20, 'sm', opts);
[V,D] = eigs(-M, 10, 0.5);
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
