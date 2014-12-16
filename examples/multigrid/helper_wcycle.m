function [u,V] = helper_wcycle(Mc, Lc, Ec, V, F, TMf2c, TMc2f, BAND, BOUNDARY, n1, n2, start, w)

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


if start == n_level
    V{start} = Mc{start} \ F{start};
else
    for i = start+1:1:n_level
        V{i} = zeros(size(F{i}));
    end
    i = start;
    for j = 1:1:n1
        V{i} = weighted_jacobi(Lc{i}, V{i}, F{i}, w);
        V{i} = Ec{i}*V{i};
    end
    res = F{i} - (Lc{i}*V{i});    
    F{i+1} = TMf2c{i}*res;
    for cnt = 1:2
        [~,V] = helper_wcycle(Mc, Lc, Ec, V, F, TMf2c, TMc2f, BAND, BOUNDARY, n1, n2, start+1, w);
    end
    v = TMc2f{i}*V{i+1};
    V{i} = V{i} + v;
    for j = 1:n2
        V{i} = jacobi(Lc{i}, V{i}, F{i});
        V{i} = Ec{i}*V{i};
    end
end

u = V{start};

end