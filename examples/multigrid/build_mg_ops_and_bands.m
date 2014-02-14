function [Mc, Lc, Ec, Eic, Rc, a_iband, a_oband, a_xcp, a_ycp, a_xg, a_yg, a_innerInOuter]= ...
    build_mg_ops_and_bands(a_band, a_xcp, a_ycp, a_xg, a_yg, a_x1d, a_y1d, p, order)

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
    [Lc{i}, Ec{i}, Rc{i}, iband, oband, a_iband{i}, a_oband{i}, a_xcp{i}, a_ycp{i}, a_xg{i}, a_yg{i}, innerInOuter] = ...
        ops_and_bands2d_test(a_x1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, a_xg{i}, a_yg{i}, a_band{i}, p, order);
    Eic{i} = Ec{i}(innerInOuter,:);
	a_innerInOuter{i} = innerInOuter;
    Mc{i} = lapsharp_unordered(Lc{i}, Ec{i}, Rc{i});
    %E = Ec{i};
    %E(innerInOuter,:) = speye(size(Eic{i}));
    %Mc{i} = Lc{i}*E;
end

end
