function [v] = gauss_jacobi(L, D, E, v, f)

w = E*v;
v = v + (f - L*w)./D;

end