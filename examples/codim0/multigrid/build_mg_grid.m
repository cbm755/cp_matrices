function [a_band, a_xcp, a_ycp, a_distg, a_bdyg, a_dx, a_x1d, a_y1d, a_xg, a_yg, a_param] = ...
    build_mg_grid(x1d_coarsest, y1d_coarsest, dx_coarsest, dx, bw, cpf, has_boundary, use_ndgrid, need_param)

% setting up cp grids, invoking function 'refine_grid'

if (nargin < 7)
    use_ndgrid = 0;
end
if (nargin < 8)
    use_ndgrid = 0;
end
if (nargin < 9)
    need_param = false;
end

n_level = round(log(dx_coarsest/dx)/log(2)) + 1;

a_dx = cell(n_level,1);
a_x1d = cell(n_level,1);
a_y1d = cell(n_level,1);
a_band = cell(n_level,1);
a_xcp = cell(n_level,1);
a_ycp = cell(n_level,1);
a_distg = cell(n_level,1);
a_bdyg = cell(n_level,1);
a_xg = cell(n_level,1);
a_yg = cell(n_level,1);
a_param = cell(n_level,1);

[xx yy] = meshgrid(x1d_coarsest, y1d_coarsest);
if need_param
    if has_boundary
        [cpx, cpy, dist, bdy, param] = cpf(xx,yy);
    else
        [cpx, cpy, dist, param] = cpf(xx,yy);
        bdy = [];
    end
else
    param = [];
    if has_boundary
        [cpx, cpy, dist, bdy] = cpf(xx,yy);
    else
        [cpx, cpy, dist] = cpf(xx,yy);
        bdy = [];
    end
end

band = find(abs(dist) <= bw*dx_coarsest);

xg = xx(band); yg = yy(band);
cpxg = cpx(band); cpyg = cpy(band);
distg = dist(band); 
if isempty(bdy)
    bdyg = []; 
else
    bdyg = bdy(band);
end
if isempty(param)
    paramg = [];
else
    paramg = param(band);
end

i = n_level;
a_dx{i} = dx_coarsest;
a_x1d{i} = x1d_coarsest;
a_y1d{i} = y1d_coarsest;
a_band{i} = band;
a_xcp{i} = cpxg;
a_ycp{i} = cpyg;
a_distg{i} = distg;
a_bdyg{i} = bdyg;
a_xg{i} = xg;
a_yg{i} = yg;
a_param{i} = paramg;

if need_param == true
for i = n_level-1:-1:1
    [a_band{i}, a_xg{i}, a_yg{i}, a_xcp{i}, a_ycp{i}, a_distg{i}, a_bdyg{i}, a_dx{i}, a_x1d{i}, a_y1d{i}] = ...
        refine_grid(1, cpf, a_dx{i+1}, a_x1d{i+1}, a_y1d{i+1}, bw, a_band{i+1}, a_distg{i+1}, a_bdyg{i+1}, use_ndgrid, need_param);
end
else
for i = n_level-1:-1:1
    [a_band{i}, a_xg{i}, a_yg{i}, a_xcp{i}, a_ycp{i}, a_distg{i}, a_bdyg{i}, a_dx{i}, a_x1d{i}, a_y1d{i}] = ...
        refine_grid(1, cpf, a_dx{i+1}, a_x1d{i+1}, a_y1d{i+1}, bw, a_band{i+1}, a_distg{i+1}, a_bdyg{i+1});
    %bw = bw*2;
end
end

end
