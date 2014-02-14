function [v] = jacobi(M, v, f)

D = diag(M);
v = v + (f - M*v) ./ D;

end