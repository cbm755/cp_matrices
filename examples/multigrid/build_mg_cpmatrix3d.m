function [Mc, Lc, Ec, L1] = build_mg_cpmatrix3d ...
         (a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order, time_dependent, dt, lambda)

% build up cp matrices M, L, E at each level of v-cycle.

n_level = length(a_band);

Ec = cell(n_level,1);
Lc = cell(n_level,1);
Mc = cell(n_level,1);

d = 3;

if nargin == 9
    
    for i = 1:1:n_level
        dx = a_x1d{i}(2) - a_x1d{i}(1);
        Ec{i} = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, p);
        Ec{i} = Ec{i}(:, a_band{i});
        Lc{i} = laplacian_3d_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i});
        Mc{i} = lapsharp(Lc{i}, Ec{i}, 2*d/dx^2);
    end
    L1 = Lc{1};
    
elseif nargin >= 10
    if nargin == 10 && time_dependent == true
        dt = [];
        lambda = 1;
    end
    if ~isempty(dt) && ~isempty(lambda)
        for i = 1:1:n_level
            Ec{i} = interp3_matrix(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, p);
            Ec{i} = Ec{i}(:, a_band{i});
            Lc{i} = laplacian_3d_matrix_test(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i});
            Mc{i} = lapsharp(Lc{i}, Ec{i});
        end
        L1 = Lc{1};
        for i = 1:1:n_level
            Mc{i} = speye(size(Mc{i})) - lambda*dt*Mc{i};
            Lc{i} = speye(size(Lc{i})) - lambda*dt*Lc{i};
        end
    end


end