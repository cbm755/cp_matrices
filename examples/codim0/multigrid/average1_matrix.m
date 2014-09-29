function [A] = average1_matrix(xx,xg,band)
        dx = xx(2) - xx(1);
        dim = 1;
        p = 2;
        N = 3;
        EXTSTENSZ = N^dim;

        Ei = repmat((1:length(xg))',1,EXTSTENSZ);
        Ej = zeros(size(Ei));

        relpt = xx(1);
        [Ibpt, Xgrid] = findGridInterpBasePt_vec(xg, p, relpt, dx);

%         % adhoc stuff for the particular test example 'poisson_1d_test',
%         % making the interpoaltion stencil total inside the interval
%         Ibpt = Ibpt - 1;
%         Xgrid = Xgrid - dx;

        weights = repmat([1/4,1/2,1/4], length(xg), 1);

        for cnt = 1:N
            Ej(:,cnt) = Ibpt + cnt - 1;
        end
 
        A = sparse(Ei(:), Ej(:), weights(:), length(xg), length(xx));
        A = A(:,band);
    end