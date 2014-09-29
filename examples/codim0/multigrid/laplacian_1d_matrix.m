 function [L] = laplacian_1d_matrix(xx,band,order)
        dx = xx(2) - xx(1);
        n = length(xx);
        e = ones(n,1);
        if order == 2
            L = spdiags([e -2*e e], -1:1, n,n);
            L = L(band,band) / dx^2;
        elseif order == 4
            L = spdiags([-e/12 4/3*e -2.5*e 4/3*e -e/12], -2:2, n,n);
            L = L(band,band) / dx^2;
        end
end