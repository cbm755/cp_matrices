function [u] = helper_vcycle_test(L, a_Ebar, a_Edouble, a_Etriple, E_out_out, E_out_in, Ecp_Omega_S, Ecp_f2c_Omega, Ecp_f2c_S, Ecp_f2c_Omega_S, V, F, FonS, TMf2c, TMc2f, BAND, BDY, n1, n2, start, w)

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
        %V{i}(bdy) = 2*F{i}(bdy) - a_Ebar{i}*V{i};
        %V{i}(bdy) = 3*F{i}(bdy) -3*a_Ebar{i}*V{i} + a_Edouble{i}*V{i}; 
        %V{i}(bdy) = 4*F{i}(bdy) - 6*a_Ebar{i}*V{i} + 4*a_Edouble{i}*V{i} - a_Etriple{i}*V{i};
    end

    % compute the residual and restrict, only make sense for interior points.
    res = F{i} - (L{i}*V{i});
    F{i+1}(~BDY{i+1}) = TMf2c{i}(~BDY{i+1},:)*res;
    
%     FonS{i} = FonS{i} - Ecp_Omega_S{i}*V{i};
%     F{i+1}(BDY{i+1}) = Ecp_f2c_Omega_S{i} * FonS{i}; %- Ecp_f2c_Omega{i} * V{i};
%     %FonS{i} = FonS{i} - Ecp_Omega_S{i}*V{i};
%     FonS{i+1} = Ecp_f2c_S{i}*FonS{i};

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
    %V{i}(bdy) = 2*F{i}(bdy) - a_Ebar{i}*V{i};
    %V{i}(bdy) = 3*F{i}(bdy) -3*a_Ebar{i}*V{i} + a_Edouble{i}*V{i}; 
    %V{i}(bdy) = 4*F{i}(bdy) - 6*a_Ebar{i}*V{i} + 4*a_Edouble{i}*V{i} - a_Etriple{i}*V{i};
    for j = 1:n2
        V{i} = weighted_jacobi(L{i}, V{i}, F{i}, w);
        %V{i} = gauss_seidel(L{i}, V{i}, F{i});
        %V{i} = jacobi(L{i}, V{i}, F{i});

        V{i}(bdy) = E_out_out{i} \ ( F{i}(bdy)  - E_out_in{i}*V{i}(nbdy) );
        %V{i}(bdy) = 2*F{i}(bdy) - a_Ebar{i}*V{i};
        %V{i}(bdy) = 3*F{i}(bdy) -3*a_Ebar{i}*V{i} + a_Edouble{i}*V{i};
        %V{i}(bdy) = 4*F{i}(bdy) - 6*a_Ebar{i}*V{i} + 4*a_Edouble{i}*V{i} - a_Etriple{i}*V{i};
    end

end

u = V{start};

end
