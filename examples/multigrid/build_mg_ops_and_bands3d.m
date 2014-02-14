function [Mc, Lc, Ec, Eic, Rc, a_iband, a_oband, a_xcp, a_ycp, a_zcp, a_xg, a_yg, a_zg, a_innerInOuter]= ...
    build_mg_ops_and_bands3d(a_band, a_xcp, a_ycp, a_zcp, a_xg, a_yg, a_zg, a_x1d, a_y1d, a_z1d, p, order)

% build up cp matrices M, L, E at each level of v-cycle.

n_level = length(a_band);

Mc = cell(n_level,1);
Lc = cell(n_level,1);
Ec = cell(n_level,1);
Eic = cell(n_level,1);
Rc = cell(n_level,1);

a_iband = cell(n_level,1);
a_oband = cell(n_level,1);
a_innerInOuter = cell(n_level,1);

for i = 1:1:n_level
    [Lc{i}, Ec{i}, Rc{i}, iband, oband, a_iband{i}, a_oband{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, a_xg{i}, a_yg{i}, a_zg{i}, innerInOuter] = ...
        ops_and_bands3d_test(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, a_xg{i}, a_yg{i}, a_zg{i}, a_band{i}, p, order);
    Eic{i} = Ec{i}(innerInOuter,:);
	a_innerInOuter{i} = innerInOuter;
    Mc{i} = lapsharp_unordered(Lc{i}, Ec{i}, Rc{i});
end

end
