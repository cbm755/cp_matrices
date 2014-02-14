function [u,err_normal_v,err_normal_f] = helper_vcycle_debug(Mc, Lc, Ec, V, F, TMf2c, TMc2f, BAND, BOUNDARY, n1, n2, start, w)

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

% seems to improve the solution of semicircle or hemisphere dealing with
% Dirichlet Boundary condtion
% V{1} = Ec{1}*V{1};

N = 2*(n_level-start);
err_normal_v = zeros(N,1);
err_normal_f = zeros(N,1);
cnt = 1;

% downward branch of the V-cycle
for i = start:1:n_level-1
%    bdy = BOUNDARY{i};
%     if ~isempty(bdy)
%         lbdy = logical(bdy);
%         F{i}(lbdy) = -F{i}(lbdy);
%     end
    for j = 1:1:n1
        V{i} = jacobi(Mc{i}, V{i}, F{i});
        
        %V{i} = weighted_jacobi(Lc{i}, V{i}, F{i},w);
        %V{i} = Ec{i}*V{i};    
    end
    %err_normal_v(cnt) = norm(V{i}-Ec{i}*V{i},inf)/norm(V{i},inf);
    err_normal_v(cnt) = norm(V{i}-Ec{i}*V{i},inf);
    disp(['pre-smoothing  level:', num2str(i),' error along normal:', num2str(err_normal_v(cnt))])
    
    % need do extension to ensure res along normal direction constant ?
    % For a cosine curve cos(t), t in [1/4, 4], Dirichlet B.C. 
    % F{i} - Ec{i}*(Lc{i}*V{i}  is better; but for Neumann B.C.
    % F{i} - Lc{i}*V{i} is better...
    %res = F{i} - Ec{i}*(Lc{i}*V{i});    
    %res = F{i} - (Lc{i}*V{i});
    
    res = F{i} - Mc{i}*V{i};
    %res = F{i} - Mc{i}*(Ec{i}*V{i});
    
    F{i+1} = TMf2c{i}*res;
    
    %err_normal_f = norm(F{i}-Ec{i}*F{i},inf)/norm(F{i+1},inf);
    err_normal_f(cnt) = norm(F{i}-Ec{i}*F{i},inf);
    %disp(['level:', num2str(i),' res along normal:', num2str(err_normal_f)])

    cnt = cnt + 1;
    % since the original rhs is constant along the normal direction, we
    % hope to keep this property at each level of V-Cycle, sometimes this
    % would make multigrid converge faster, but sometimes less accurate.
    % Seems important when dealing with Dirichlet B.C.
    % F{i+1} = Ec{i+1}*F{i+1};    
end

% solve on the coarsest grid
if isa(Mc, 'cell')
    V{n_level} = Mc{n_level} \ F{n_level};
    %V{n_level} = Ec{n_level}*V{n_level};
    %[V{n_level} flag] = gmres(Mc{n_level}, F{n_level}, [], 1e-6);
else
    %V{n_level} = Mc \ F{n_level};
    [V{n_level} flag] = gmres(Mc, F{n_level}, [], 1e-6);
end
% [V{n_level} flag] = gmres(Mc{n_level}, F{n_level}, [], 1e-10);
% [V{n_level} flag] = gmres(Mc{n_level}, F{n_level}, 5, 1e-6);

% upward branch of the V-cycle
for i = n_level-1:-1:start   
    
    v = TMc2f{i}*V{i+1};
    % here Ec{i}*v will improve the accuracy when using the Crank-Nicolson
    % Scheme; and Ec{i}*v is essentially neccessary in solving poisson
    % euqation on a semicircle with Neumann and Dirichlet Boundary conditions; 
    % and Ec{i}*v is also important for Dirichlet B.C. on hemisphere
    %V{i} = V{i} + Ec{i}*v;
    V{i} = V{i} + v;
    for j = 1:n2
        % we should use Ec{i}*V{i} as the initial value
        V{i} = jacobi(Mc{i}, V{i}, F{i});
        
        %V{i} = weighted_jacobi(Lc{i}, V{i}, F{i}, w);
        %V{i} = Ec{i}*V{i};
        
    end
    %err_normal_v = norm(V{i}-Ec{i}*V{i},inf)/norm(V{i},inf);
    err_normal_v(cnt) = norm(V{i}-Ec{i}*V{i},inf);
    disp(['post-smoothing level:', num2str(i),' error along normal:', num2str(err_normal_v(cnt))])

    res = F{i} - Mc{i}*V{i};
    %res = F{i} - Mc{i}*(Ec{i}*V{i});
    err_normal_f(cnt) = norm(res-Ec{i}*res,inf);
    
    cnt = cnt + 1;
end

u = V{start};

end
