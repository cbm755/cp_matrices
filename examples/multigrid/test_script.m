
%% trying to find whether there are better lambda's to make those rows having non-negative positive entries.
[i,j] = ind2sub(size(flag2),find(flag2));
ii = unique(i);

cnt = 0;
sad = 0;
vec = zeros(size(ii));
for k = 1:length(ii)
    flag = false;
    for lambda = 2:0.1:6
        e_k = zeros(1,length(E));
        e_k(ii(k)) = 1;
        tmp = E(ii(k),:)*Lc{n_level} - lambda/dx_coarsest^2*(e_k-Ec{n_level}(ii(k),:));
        
%         ind = (tmp<0);
%         tmp(ind)
        
        if (nnz(tmp<-1e-10) == 1)
%              [k, lambda]
            cnt = cnt + 1;
            flag = true;
            break
        end    
    end
    if ~flag
        k
        sad = sad+1;
        vec(sad) = k;
    end
end
vec = vec(1:sad);
% what are these hopeless points?
[a_xcp{1}(ii(vec)) a_ycp{1}(ii(vec)) a_zcp{1}(ii(vec))]


%% following codes aim to look at where are these hopeless points on the sphere.
p = 1;
ddx = [dx_coarsest dx_coarsest dx_coarsest];
[Ibpt,Xgrid] = findGridInterpBasePt_vec({a_xcp{1}(ii(vec)) a_ycp{1}(ii(vec)) a_zcp{1}(ii(vec))},p,[a_x1d{1},a_y1d{1},a_z1d{1}],ddx);

plot3(Xgrid{1},Xgrid{2},Xgrid{3},'.')
hold on, plot3(a_xcp{1}(ii(vec)),a_ycp{1}(ii(vec)),a_zcp{1}(ii(vec)),'r.')
xlabel('x'),ylabel('y'),zlabel('z')
grid on, axis equal


diagGamma = diag(GAMMA);
coeff = diagGamma(ii)*0.2^2;

k = 2;
[a_xcp{1}(ii(k)), a_ycp{1}(ii(k)), a_zcp{1}(ii(k))]
[Xgrid{1}(k),Xgrid{2}(k),Xgrid{3}(k)]



