

dx = 0.1;                   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;
w1d = x1d;

nx = length(x1d);
ny = length(y1d);
nz = length(z1d);
nw = length(w1d);

[xx yy zz ww] = ndgrid(x1d, y1d, z1d, w1d);
[th, r] = cart2pol(xx, yy);
[cpx, cpy] = pol2cart(th, 1);
cpz = zeros(size(cpx));
cpw = zeros(size(cpx));
dist = sqrt((xx-cpx).^2 + (yy-cpy).^2 + (zz-cpx).^2 + (ww-cpw).^2);


%% Banding: do calculation in a narrow band around the sphere
dim = 4;  % dimension
p = 3;    % interpolation order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpx = cpx(band); cpy = cpy(band); cpz = cpz(band); cpw = cpw(band);
xg = xx(band); yg = yy(band); zg = zz(band); wg = ww(band);


%E1 = interp3_matrix(x1d,y1d,z1d,  cpx, cpy, cpz, 3, band);
E = interpn_matrix({x1d y1d z1d w1d}, [cpx cpy cpz cpw], 3, band);