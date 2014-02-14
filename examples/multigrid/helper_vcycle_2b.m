function [u] = helper_vcycle_2b(Mc, Lc, Ec, Eic, D, V, F, TMf2c, TMc2f, BAND, BOUNDARY, n1, n2, start, w)

% V: a vector storing solution vectors of each V-cycle level
% M: a matrix storing iterative matrices of each V-cycle level
% F: a vector storing the right hand side (residual) of each V-cycle level
% TMf2c:   transform matrices from fine grid to coarse grid
% TMc2f:   transform matrices from coarse grid to fine grid
% p:  degree of interpolation polynomial
% n1: number of gauss-seidel iterations of the downward branch of the V-cycle
% n2: number of gauss-seidel iterations of the upward branch of the V-cycle

% number of the levels of multigrid
n_level = length(BAND);


% downward branch of the V-cycle
for i = start:1:n_level-1
    
    for j = 1:1:n1
    V{i} = gauss_jacobi(Lc{i}, D{i}, Ec{i}, V{i}, F{i});
    %V{i} = gauss_seidel(Mc{i}, V{i}, F{i});
    %V{i} = jacobi(Mc{i}, V{i}, F{i});
    V{i} = Eic{i}*V{i};
    end
    % the following step of extension seems very import to make soultion
    % accurate, but not know why
    % V{i} = Eic{i}*V{i};
    
    %res = F{i} - Mc{i}*V{i};
    
    res = F{i} - Lc{i}*(Ec{i}*V{i});
    F{i+1} = TMf2c{i}*res;
    
end

% solve on the coarsest grid
 V{n_level} = Mc{n_level} \ F{n_level};
% [V{n_level} flag] = gmres(Mc{n_level}, F{n_level}, [], 1e-10);
% [V{n_level} flag] = gmres(Mc{n_level}, F{n_level});

% upward branch of the V-cycle
for i = n_level-1:-1:start   
    
    v = TMc2f{i}*V{i+1};
    V{i} = V{i} + v;
    
    for j = 1:n2
        V{i} = gauss_jacobi(Lc{i}, D{i}, Ec{i}, V{i}, F{i});
        %V{i} = gauss_seidel(Mc{i}, V{i}, F{i});
        %V{i} = jacobi(Mc{i}, V{i}, F{i});
        V{i} = Eic{i}*V{i};
    end
    %V{i} = Eic{i}*V{i};
end

u = V{start};

end
