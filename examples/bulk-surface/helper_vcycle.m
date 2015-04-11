function [u, v] = helper_vcycle(Au, Lu, Eu, A_coarsest, Av, Eoo_v, Eoi_v, Bvu, Buv, U, V, G, F, TMf2c, TMc2f, TMf2c_S, TMc2f_S, ...
            a_band, a_band_S, a_bdyg, n1, n2, start, w)

% number of the levels of multigrid
n_level = length(a_band);

for i = start+1:1:n_level
    F{i}(a_bdyg{i}) = 0;
end

% downward branch of the V-cycle
for i = start:1:n_level-1
    bdy = a_bdyg{i};
    nbdy = ~bdy;
    for j = 1:1:n1
        %for cnt = 1:n1
        % Fix u, do relaxation and modify ghost points values for v
        V{i} = weighted_jacobi(Av{i}, V{i}, F{i}, w);
        V{i}(bdy) = Eoo_v{i} \ ( -Bvu{i}*U{i}  - Eoi_v{i}*V{i}(nbdy) );
        %end
        
        %for cnt = 1:n1
        % Fix v, do relaxation and closest point extension for u
        U{i} = weighted_jacobi(Au{i}, U{i}, G{i}-Buv{i}*V{i}, 1);
        U{i} = Eu{i}*U{i};
        %end
    end

    % compute the residual for v and restrict, only make sense for interior points.
    res_v = F{i} - (Av{i}*V{i});
    F{i+1}(~a_bdyg{i+1}) = TMf2c{i}(~a_bdyg{i+1},:)*res_v;
    % compute residual for u and restrict
    res_u = G{i} - Buv{i}*V{i} - (Au{i}*(U{i}));
    G{i+1} = TMf2c_S{i}*res_u;
end

% solve on the coarsest grid
soln = A_coarsest \ [F{n_level}; G{n_level}];
V{n_level} = soln(1:length(a_band{n_level}));
U{n_level} = soln(length(a_band{n_level})+1:end);

% upward branch of the V-cycle
for i = n_level-1:-1:start   
    bdy = a_bdyg{i}; 
    nbdy = ~bdy;

    u = TMc2f_S{i}*U{i+1};
    U{i} = U{i} + Eu{i}*u;
    
    v = TMc2f{i}*V{i+1};
    V{i} = V{i} + v;
    V{i}(bdy) = Eoo_v{i} \ ( -Bvu{i}*U{i}  - Eoi_v{i}*V{i}(nbdy) );

    for j = 1:n2
        %for cnt = 1:n2
        % Fix u, do relaxation and modify ghost points values for v
        V{i} = weighted_jacobi(Av{i}, V{i}, F{i}, w);
        V{i}(bdy) = Eoo_v{i} \ ( -Bvu{i}*U{i}  - Eoi_v{i}*V{i}(nbdy) );
        %end
        
        %for cnt = 1:n2
        % Fix v, do relaxation and closest point extension for u
        U{i} = weighted_jacobi(Au{i}, U{i}, G{i}-Buv{i}*V{i}, 1);
        U{i} = Eu{i}*U{i};
        %end
    end

end

u = U{start};
v = V{start};
end
