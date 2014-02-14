function [TMf2c, TMc2f] = helper_set_TM(X, Y, XCP, YCP, BAND, BDYG, p_f2c, p_c2f)

%% construct the transfrom matrices to do the restriction and prolongation between the neibouring two levels of the V-cycle

% TMf2c:   a vector storing transform matrices from fine grid to coarse grid, 
%          each matrix multiplying the fine grid's vector will restrict the vector to the coarse grid

% TMc2f:   similar to TMf2c, doing the prolongation from coarse grid to fine grid

% X,Y:     x & y coordinates which define the whole 2D domain at each level of the V-cycle
% XCP,YCP: coordinates of all closest points of each level of the V-cycle
% BAND:    the innerbands of all levels of the V-cycle


if (nargin < 10)
    use_ndgrid = false;
end

n_level = length(BAND);

TMf2c = cell(n_level,1);
TMc2f = cell(n_level,1);

for i = 1:1:n_level-1
    x = X{i};
    y = Y{i};
    xi = XCP{i+1};
    yi = YCP{i+1};
    band = BAND{i};
    bdyg = BDYG{i+1};
    T = interp2_matrix( x, y, xi, yi, p_f2c);
    T = T(:,band);
    bdyg = logical(bdyg);
    %T(bdyg,:) = -T(bdyg,:);

    %TMf2c{i} = T*A{i};
    TMf2c{i} = T;
end


for i = 1:1:n_level-1
    x = X{i+1};
    y = Y{i+1};
    xi = XCP{i};
    yi = YCP{i};
    band = BAND{i+1};
    bdyg = BDYG{i};

    T = interp2_matrix( x, y, xi, yi, p_c2f);
    T = T(:,band);
    bdyg = logical(bdyg);
    %T(bdyg,:) = -T(bdyg,:);
    TMc2f{i} = T;

end

end
