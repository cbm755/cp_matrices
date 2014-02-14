function [Lc, Ec, L1, M1, Mn] = build_mg_cpmatrix3d_test ...
         (a_band, a_xcp, a_ycp, a_zcp, a_x1d, a_y1d, a_z1d, p, order, time_dependent, dt, lambda)

% build up cp matrices M, L, E at each level of v-cycle.

n_level = length(a_band);

Ec = cell(n_level,1);
Lc = cell(n_level,1);

if nargin == 9
    
    for i = 1:1:n_level
        Ec{i} = interp3_matrix_test(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, p);
        Ec{i} = Ec{i}(:, a_band{i});
        Lc{i} = laplacian_3d_matrix_test(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i});
    end
    L1 = Lc{1};
    M1 = lapsharp(Lc{1},Ec{1});
    Mn = lapsharp(Lc{n_level},Ec{n_level});
    
elseif nargin >= 10
    if nargin == 10 && time_dependent == true
        dt = [];
        lambda = 1;
    end
    if ~isempty(dt) && ~isempty(lambda)
        for i = 1:1:n_level
            Ec{i} = interp3_matrix_test(a_x1d{i}, a_y1d{i}, a_z1d{i}, a_xcp{i}, a_ycp{i}, a_zcp{i}, p);
            Ec{i} = Ec{i}(:, a_band{i}); 
        end

        L1 = laplacian_3d_matrix_test(a_x1d{1}, a_y1d{1}, a_z1d{1}, order, a_band{1}, a_band{1});
        tic; M1 = lapsharp(L1, Ec{1}); toc;

        for i = 1:1:n_level-1
            Lc{i} = laplacian_3d_matrix_test(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i}, ...
                                             time_dependent, dt, lambda);
        end
        i = n_level;
        Lc{i} = laplacian_3d_matrix_test(a_x1d{i}, a_y1d{i}, a_z1d{i}, order, a_band{i}, a_band{i});
        tic; Mn = lapsharp(Lc{i},Ec{i}); toc;
    end
end


end