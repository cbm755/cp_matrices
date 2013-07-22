function [pass, str] = test_vectorized_interp3_matrix()
  str = 'test vectorized interp3 method';

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
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
xg = xx(band); yg = yy(band); zg = zz(band);




%% Construct an interpolation matrix for closest point
disp('Constructing interpolation matrix');
T=cputime();
E2full = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, p);
E2 = E2full(:,band);
cputime()-T

disp('Constructing interpolation matrix w/ banding');
T=cputime();
E3 = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, p, band);
cputime()-T

c = c + 1;
pass(c) = max(max(abs(E2-E3))) == 0;


disp('Constructing interpolation matrix, index form');
%profile on
T=cputime();
[E4i, E4j, E4s] = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, p);
cputime()-T
%profile off
%profile viewer


disp('Constructing interpolation matrix with loops');
tic
Efull = interp3_matrix_oldloop(x1d,y1d,z1d, cpxg,cpyg,cpzg, p);
E = Efull(:,band);
toc

c = c + 1;
pass(c) = max(max(abs(E-E2))) == 0;

%% compare to the extracted contents of E
% no guarantee order doesn't change, and there is the "geometric vs
% linalg" column-rank of E issue.
[I,J,V] = find(E);
f = E4s == 0;
c = c + 1;
pass(c) = max(max(abs(sort(E4s(~f)) - sort(V)))) == 0;
