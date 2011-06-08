%% Eigenvalue problem on hemisphere
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

% This example computes eigenvalues and eigenfunctions of the
% Laplace--Beltrami operator on a hemisphere.  It also demonstates how
% to impose 2nd-order accurate Neumann and Dirichlet BCs [Macdonald,
% Brandman, Ruuth 2011].


% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');

dx = 0.2/2^1;  % grid spacing
R=1; %Radius of the circle

pad = 5;
x1d=(-R-pad*dx):dx:(R+pad*dx)';
y1d=(-R-pad*dx):dx:(R+pad*dx)';
z1d=(-R-pad*dx):dx:(R+pad*dx)';
nx=length(x1d);
ny=length(y1d);
nz=length(z1d);
[x3d,y3d,z3d]=meshgrid(x1d,y1d,z1d);

% using the standard CP function, we get a homogeneous Neuamann BC
% [Ruuth & Merriman 2008]
%[cpx,cpy,cpz, dist, bdy] = cpHemisphere(x3d,y3d,z3d R);

% Using "cpbar" [Macdonald, Brandman, Ruuth 2011]:
cpf = @cpHemisphere;
[cpx,cpy,cpz, dist, bdy] = cpbar_3d(x3d,y3d,z3d, cpf, R);

dim = 3;
p = 3;  % degree interp
order = 2;  % laplacian order
bw = 1.001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
if (pad < bw)
  warning('maybe insufficient padding');
end
band = find(dist <= bw*dx);
outband = find(dist > bw*dx);

% Last place we need the 3D array
xg = x3d(band);
yg = y3d(band);
zg = z3d(band);
cpxg = cpx(band);
cpyg = cpy(band);
cpzg = cpz(band);
bdyg = bdy(band);


%% discrete operators
disp('building laplacian and interp matrices');
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band,band);
E = interp3_matrix_band(x1d,y1d,z1d, cpxg, cpyg, cpzg, p, band);
% dirichlet BCs
%E(bdyb,:) = -E(bdyb,:);
% iCPM matrix
M = diagSplit(L,E);


%% plotting grid
[xp,yp,zp] = paramHemisphere(48, R);
%[xp,yp,zp] = paramSphere(32, R);
xp1 = xp(:);  yp1 = yp(:);  zp1 = zp(:);
disp('building plotting matrices');
Eplot = interp3_matrix_band(x1d,y1d,z1d, xp1,yp1,zp1, p, band);
Eplot0 = interp3_matrix_band(x1d,y1d,z1d, xp1,yp1,zp1, 0, band);




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


tic
figure(2); clf;
for i=1:4
  lambda = Lambda(i);
  eigenvec = V(:,I(i));

  uplot = Eplot*E*real(eigenvec);
  uplot = reshape(uplot, size(xp));
  subplot(2,2,i);
  surf(xp,yp,zp, uplot);

  xlabel('x'); ylabel('y'); zlabel('z');
  title(['eigenvalue = ' num2str(lambda)]);
  axis equal
  shading interp
  camlight left
end
plottime = toc