function [Mc, Lc, Ec] = build_mg_cpmatrix(a_band, a_xcp, a_ycp, a_x1d, a_y1d, p, order)

% build up cp matrices M, L, E at each level of v-cycle.

n_level = length(a_band);

Ec = cell(n_level,1);
Lc = cell(n_level,1);
Mc = cell(n_level,1);

d = 2;

for i = 1:1:n_level
    dx = a_x1d{i}(2) - a_x1d{i}(1);
    Ec{i} = interp2_matrix(a_x1d{i}, a_y1d{i}, a_xcp{i}, a_ycp{i}, p);
    Ec{i} = Ec{i}(:, a_band{i});
    Lc{i} = laplacian_2d_matrix(a_x1d{i}, a_y1d{i}, order, a_band{i}, a_band{i});
    %Mc{i} = Lc{i}*Ec{i} - 2*d/dx^2*(speye(size(Ec{i}))-Ec{i});
    Mc{i} = lapsharp(Lc{i}, Ec{i}, 2*d/dx^2);
    %Mc{i} = lapsharp(Lc{i}, Ec{i});
end

%Mc{n_level} = lapsharp(Lc{i}, Ec{i}, 2*d/dx);

end
