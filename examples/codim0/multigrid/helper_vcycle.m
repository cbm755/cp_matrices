function [u] = helper_vcycle(L, E_out_out, E_out_in, V, F, TMf2c, TMc2f, BAND, BDY, n1, n2, start, w)

% number of the levels of multigrid
n_level = length(BAND);

for i = start+1:1:n_level
    F{i}(BDY{i}) = 0;
end

% downward branch of the V-cycle
for i = start:1:n_level-1
    bdy = BDY{i};
    nbdy = ~bdy;
    for j = 1:1:n1
        % do one step of Jacobi iteration, only make sense for interior points.
        V{i} = weighted_jacobi(L{i}, V{i}, F{i}, w);
        %V{i} = gauss_seidel(L{i}, V{i}, F{i});
        %V{i} = jacobi(L{i}, V{i}, F{i});

        % modify values of ghost points.
        V{i}(bdy) = E_out_out{i} \ ( F{i}(bdy)  - E_out_in{i}*V{i}(nbdy) );
    end

    % compute the residual and restrict, only make sense for interior points.
    res = F{i} - (L{i}*V{i});
    F{i+1}(~BDY{i+1}) = TMf2c{i}(~BDY{i+1},:)*res;

end

% solve on the coarsest grid
V{n_level} = L{n_level} \ F{n_level};
% [V{n_level} flag] = gmres(L{n_level}, F{n_level}, [], 1e-10);
% [V{n_level} flag] = gmres(L{n_level}, F{n_level}, 5, 1e-6);

% upward branch of the V-cycle
for i = n_level-1:-1:start   
    bdy = BDY{i}; 
    nbdy = ~bdy;
    
    v = TMc2f{i}*V{i+1};
    V{i} = V{i} + v;
    % modify values of ghost points.
    V{i}(bdy) = E_out_out{i} \ ( F{i}(bdy)  - E_out_in{i}*V{i}(nbdy) );

    for j = 1:n2
        V{i} = weighted_jacobi(L{i}, V{i}, F{i}, w);
        %V{i} = gauss_seidel(L{i}, V{i}, F{i});
        %V{i} = jacobi(L{i}, V{i}, F{i});

        V{i}(bdy) = E_out_out{i} \ ( F{i}(bdy)  - E_out_in{i}*V{i}(nbdy) );

    end

end

u = V{start};

end
