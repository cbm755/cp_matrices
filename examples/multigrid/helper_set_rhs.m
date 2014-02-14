function [V, F] = helper_set_rhs(a_xcp, a_ycp, rhsfn, start)

n_level = length(a_xcp);
F = cell(n_level,1);
V = cell(n_level,1);

for i = start:1:n_level
    [thg, rg] = cart2pol(a_xcp{i}, a_ycp{i});
    F{i} = rhsfn(thg, rg);
    %F{i} = zeros(size(thg));
    
    V{i} = zeros(size(F{i}));    
    %V{i} = 0.5*ones(size(F{i}));
    %V{i} = rand(size(F{i}));
end