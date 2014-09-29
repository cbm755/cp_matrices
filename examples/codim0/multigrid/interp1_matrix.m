function [E] = interp1_matrix(xx,xg,p,band)
        dx = xx(2) - xx(1);
        dim = 1;
        N = p+1;
        EXTSTENSZ = N^dim;

        Ei = repmat((1:length(xg))',1,EXTSTENSZ);
        Ej = zeros(size(Ei));

        relpt = xx(1);
        [Ibpt, Xgrid] = findGridInterpBasePt_vec(xg, p, relpt, dx);

%         % adhoc stuff for the particular test example 'poisson_1d_test',
%         % making the interpoaltion stencil total inside the interval
%         Ibpt = Ibpt - 1;
%         Xgrid = Xgrid - dx;

        weights = LagrangeWeights1D_vec(Xgrid, xg, dx, N);

        for cnt = 1:N
            Ej(:,cnt) = Ibpt + cnt - 1;
        end
 
        E = sparse(Ei(:), Ej(:), weights(:), length(xg), length(xx));
        E = E(:,band);
    end