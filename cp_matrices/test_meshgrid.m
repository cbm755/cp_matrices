function [weights, weights1] = test_meshgrid(N, X, p, relpt, dx)
% ii = I(1) + (0:N-1);
% jj = I(2) + (0:N-1);
% [iig, jjg] = meshgrid(ii,jj);
% iig = iig(:);
% jjg = jjg(:);
% 
% ii1 = (0:N-1);
% jj1 = (0:N-1);
% [iig1, jjg1] = meshgrid(ii1,jj1);
% iig1 = iig1(:) + I(1);
% jjg1 = jjg1(:) + I(2);
dim = length(dx);
[I, Xgrid] = findGridInterpBasePt(X, p, relpt, dx);

if dim == 2
    xweights = LagrangeWeights1D(Xgrid(1), X(1), dx(1), N);
    yweights = LagrangeWeights1D(Xgrid(2), X(2), dx(2), N);
    [xweights_tmp, yweights_tmp] = meshgrid(xweights, yweights);
    weights = xweights_tmp .* yweights_tmp;
    weights = weights(:);

    weights1 = zeros(1,N^2);
    for i = 1:N
        for j = 1:N
            weights1(:,N*(i-1)+j) = xweights(:,i) .* yweights(:,j);
        end
    end

elseif dim == 3
    xw = LagrangeWeights1D(Xgrid(1), X(1), dx(1), N);
    yw = LagrangeWeights1D(Xgrid(2), X(2), dx(2), N);
    zw = LagrangeWeights1D(Xgrid(3), X(3), dx(3), N);
    [xw_tmp, yw_tmp, zw_tmp] = meshgrid(xw, yw, zw);
    weights = xw_tmp .* yw_tmp .* zw_tmp;
    weights = weights(:);
   
    weights1 = zeros(1,N^3);
    for i = 1:N
        for j = 1:N
            for k = 1:N
                weights1(:,N^2*(k-1)+N*(i-1)+j) = xw(:,i) .* yw(:,j) .* zw(:,k);
            end
        end
    end
    
end

end
