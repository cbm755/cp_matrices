function [pass, str] = test_vectorized_fd3d()
  str = 'test vectorized finite differences in 3d';

  pass = [];
  c = 0;

dx = 0.1;
% make vectors of x, y, positions of the grid
x1d = ((-1-5*dx):dx:(1+5*dx))';
y1d = ((-1-6*dx):dx:(1+6*dx))';
z1d = ((-1-7*dx):dx:(1+7*dx))';
nx = length(x1d);
ny = length(y1d);
nz = length(z1d);



%% Find closest points on the surface
[xx yy zz] = meshgrid(x1d, y1d, z1d);
[cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);


%% Banding based on Ruuth--Merriman distance
dim = 3;  % dimension
p = 3;    % interpolation order
order = 2;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
xg = xx(band); yg = yy(band); zg = zz(band);




%% Construct an interpolation matrix for closest point
disp('Constructing interpolation matrix');
T=tic;
%profile on
L2 = laplacian_3d_matrix(x1d, y1d, z1d, order, band);
%profile off
%profile viewer
%E2 = E2full(:,band);
toc(T)

%disp('Constructing interpolation matrix w/ banding');
%T=tic;
%L2E3 = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, p, band);
%[E4i, E4j, E4s] = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, p);
%toc(T)

%c = c + 1;
%pass(c) = max(max(abs(E2-E3))) == 0;

disp('Constructing laplacian matrix with yujia''s repmat');
T=tic;
L3 = laplacian_3d_matrix_test(x1d, y1d, z1d, order, band);
toc(T)



disp('Constructing laplacian matrix with loops');
T=tic;
L = laplacian_3d_matrix(x1d, y1d, z1d, order, band, band, 0, 1);
toc(T)


c = c + 1;
pass(c) = max(max(abs(L-L2))) == 0;

c = c + 1;
pass(c) = max(max(abs(L2-L3))) == 0;
