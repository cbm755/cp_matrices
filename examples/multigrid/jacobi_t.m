function [v] = jacobi_t(M, DM, v, f)

v = v + (f - M*v) ./ DM;

end