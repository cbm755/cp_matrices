function [a_band, a_xcp, a_distg, a_bdyg, a_dx, a_x1d, a_xg] = ...
    build_mg_grid_1d(x1d_coarsest, dx_coarsest, dx, bw, cpf)

n_level = round(log(dx_coarsest/dx)/log(2)) + 1;

a_dx = cell(n_level,1);
a_x1d = cell(n_level,1);
a_band = cell(n_level,1);
a_xcp = cell(n_level,1);
a_distg = cell(n_level,1);
a_bdyg = cell(n_level,1);
a_xg = cell(n_level,1);
        
left = x1d_coarsest(1);
right = x1d_coarsest(end);
for i = n_level:-1:1
    a_dx{i} = dx_coarsest / 2^(n_level-i);
    a_x1d{i} = left:a_dx{i}:right;
    [cpx, dist, bdy] = cpf(a_x1d{i});
    band = find(abs(dist) <= bw*a_dx{i});

    a_band{i} = band;
    a_xg{i} = a_x1d{i}(band); 
    a_xcp{i} = cpx(band);
    a_distg{i} = dist(band); 
    a_bdyg{i} = bdy(band);
end

end
