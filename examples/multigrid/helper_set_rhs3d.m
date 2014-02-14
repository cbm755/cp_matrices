function [V, F] = helper_set_rhs3d(a_xcp, a_ycp, a_zcp, rhsfn, start)

n_level = length(a_xcp);
F = cell(n_level,1);
V = cell(n_level,1);

for i = start:1:n_level
    [th, phi, R] = cart2sph(a_xcp{i}, a_ycp{i}, a_zcp{i});
    F{i} = rhsfn(th, phi, R);
    V{i} = zeros(size(F{i}));
end