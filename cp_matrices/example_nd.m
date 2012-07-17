

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

%[xx yy zz] = meshgrid(x1d, y1d, z1d);
[xx yy zz] = ndgrid(x1d, y1d, z1d);
[th, r] = cart2pol(xx, yy);
[cpx, cpy] = pol2cart(th, 1);
cpz = zeros(size(cpx));
dist = sqrt((xx-cpx).^2 + (yy-cpy).^2 + (zz-cpx).^2);
%cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);


%% Banding: do calculation in a narrow band around the sphere
dim = 3;  % dimension
p = 3;    % interpolation order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
xg = xx(band); yg = yy(band); zg = zz(band);


E1 = interp3_matrix(x1d,y1d,z1d,  cpx, cpy, cpz, 3, band, true);
E2 = interpn_matrix({x1d y1d z1d}, [cpx cpy cpz], 3, band);

%
E1-E2

[Dxb,Dxf, Dyb,Dyf, Dzb,Dzf] = ...
  firstderiv_upw1_3d_matrices(x1d,y1d,z1d, band, band, true);

[Db,Df] = firstderiv_upw1_nd_matrices({x1d,y1d,z1d}, band);

Dxb - Db{1}
Dyb - Db{2}
Dzb - Db{3}

Dxf - Df{1}
Dyf - Df{2}
Dzf - Df{3}


[Dxc,Dyc,Dzc] = firstderiv_cen2_3d_matrices(x1d,y1d,z1d, band, band, true);

Dc = firstderiv_cen2_nd_matrices({x1d,y1d,z1d}, band);

Dxc - Dc{1}
Dyc - Dc{2}
Dzc - Dc{3}

[Dxxc,Dyyc,Dzzc] = secondderiv_cen2_3d_matrices(x1d,y1d,z1d, band, band, true);

D2c = secondderiv_cen2_nd_matrices({x1d,y1d,z1d}, band);


L2 = laplacian_3d_matrix(x1d,y1d,z1d, 2, band, band, true);
% TODO
%L2n = laplacian_nd_matrix({x1d,y1d,z1d}, 2, band, band, true);


L4 = laplacian_3d_matrix(x1d,y1d,z1d, 4, band, band, true);
%L4n = laplacian_nd_matrix({x1d,y1d,z1d}, 4, band, band, true);
