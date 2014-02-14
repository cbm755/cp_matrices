function [a_band, a_xcp, a_ycp, a_zcp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_z1d, a_xg, a_yg, a_zg] = ...
    build_mg_cpgrid3d(x1d_coarsest, y1d_coarsest, z1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary)

% setting up cp grids, invoking function 'refine_grid'

n_level = round(log(dx_coarsest/dx)/log(2)) + 1;

a_dx = cell(n_level,1);
a_x1d = cell(n_level,1);
a_y1d = cell(n_level,1);
a_z1d = cell(n_level,1);
a_band = cell(n_level,1);
a_xcp = cell(n_level,1);
a_ycp = cell(n_level,1);
a_zcp = cell(n_level,1);
a_distg = cell(n_level,1);
a_bdyg = cell(n_level,1);
a_xg = cell(n_level,1);
a_yg = cell(n_level,1);
a_zg = cell(n_level,1);

[xx yy zz] = meshgrid(x1d_coarsest, y1d_coarsest, z1d_coarsest);
tic
if has_boundary
    [cpx, cpy, cpz, dist, bdy] = cpf(xx(:),yy(:),zz(:));
else
    [cpx, cpy, cpz, dist] = cpf(xx(:),yy(:),zz(:));
    bdy = [];
end
coarsest_level_cp_time = toc

band = find(abs(dist) <= bw*dx_coarsest);

xg = xx(band); yg = yy(band); zg = zz(band);
cpxg = cpx(band); cpyg = cpy(band); cpzg = cpz(band);
distg = dist(band); 
if isempty(bdy)
    bdyg = []; 
else
    bdyg = bdy(band);
end

i = n_level;
a_dx{i} = dx_coarsest;
a_x1d{i} = x1d_coarsest;
a_y1d{i} = y1d_coarsest;
a_z1d{i} = z1d_coarsest;
a_band{i} = band;
a_xcp{i} = cpxg;
a_ycp{i} = cpyg;
a_zcp{i} = cpzg;
a_distg{i} = distg;
a_bdyg{i} = bdyg;
a_xg{i} = xg;
a_yg{i} = yg;
a_zg{i} = zg;

for i = n_level-1:-1:1
    [a_band{i}, a_xg{i}, a_yg{i}, a_zg{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, a_distg{i}, a_bdyg{i}, a_dx{i}, a_x1d{i}, a_y1d{i}, a_z1d{i}] = ...
        refine_grid(1, cpf, a_dx{i+1}, a_x1d{i+1}, a_y1d{i+1}, a_z1d{i+1}, bw, a_band{i+1}, a_distg{i+1}, a_bdyg{i+1});
end

end
