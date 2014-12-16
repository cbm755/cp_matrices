function [v] = weighted_jacobi(M, v, f, w)

D = diag(M);
v = v + w*(f - M*v)./D;

end