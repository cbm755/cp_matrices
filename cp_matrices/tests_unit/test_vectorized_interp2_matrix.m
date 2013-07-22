function [pass, str] = test_vectorized_interp2_matrix()
  str = 'test vectorized interp2 method';

  pass = [];
  c = 0;

dx = 0.01;
% make vectors of x, y, positions of the grid
x1d = ((-1-5*dx):dx:(1+5*dx))';
y1d = ((-1-6*dx):dx:(1+6*dx))';
nx = length(x1d);
ny = length(y1d);


%% Find closest points on the surface
[xx yy] = meshgrid(x1d, y1d);
[cpx, cpy, dist] = cpCircle(xx, yy);
cpxg = cpx(:); cpyg = cpy(:);


%% Banding based on Ruuth--Merriman distance
dim = 2;  % dimension
p = 3;    % interpolation order
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band);
xg = xx(band); yg = yy(band);



%% Construct an interpolation matrix for closest point
disp('Constructing interpolation matrix');
T=cputime();
E2full = interp2_matrix(x1d, y1d, cpxg, cpyg, p);
E2 = E2full(:,band);
cputime()-T

disp('Constructing interpolation matrix w/ banding');
T=cputime();
E3 = interp2_matrix(x1d, y1d, cpxg, cpyg, p, band);
cputime()-T

c = c + 1;
pass(c) = max(max(abs(E2-E3))) == 0;


disp('Constructing interpolation matrix, index form');
%profile on
T=cputime();
[E4i, E4j, E4s] = interp2_matrix(x1d, y1d, cpxg, cpyg, p);
cputime()-T
%profile off
%profile viewer


disp('Constructing interpolation matrix with loops');
tic
Efull = interp2_matrix_oldloop(x1d, y1d, cpxg, cpyg, p);
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
