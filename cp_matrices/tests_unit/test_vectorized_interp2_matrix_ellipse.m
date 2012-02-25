function [pass, str] = test_vectorized_interp2_matrix_ellipse()
  str = 'test vectorized interp2 method (ellipse)';

  pass = [];
  c = 0;

dx = 0.1;
% make vectors of x, y, positions of the grid
x1d = ((-1.5-5*dx):dx:(1.5+5*dx))';
y1d = ((-1-6*dx):dx:(1+6*dx))';
nx = length(x1d);
ny = length(y1d);


%% Find closest points on the surface
[xx yy] = meshgrid(x1d, y1d);
[cpx, cpy, dist] = cpEllipse(xx, yy);
cpxg = cpx(:); cpyg = cpy(:);


%% Banding based on Ruuth--Merriman distance
dim = 2;  % dimension
p = 0;    % interpolation order
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band);
xg = xx(band); yg = yy(band);


E2 = interp2_matrix(x1d, y1d, cpxg, cpyg, p, band);

Efull = interp2_matrix_oldloop(x1d, y1d, cpxg, cpyg, p);
E = Efull(:,band);

c = c + 1;
pass(c) = max(max(abs(E-E2))) == 0;

