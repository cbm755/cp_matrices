function [TMf2c, TMc2f] = helper_set_TM_3d(X, Y, Z, XG, YG, ZG, BAND, p_f2c, p_c2f)

%% construct the transfrom matrices to do the restriction and prolongation between the neibouring two levels of the V-cycle

% TMf2c:   a vector storing transform matrices from fine grid to coarse grid, 
%          each matrix multiplying the fine grid's vector will restrict the vector to the coarse grid

% TMc2f:   similar to TMf2c, doing the prolongation from coarse grid to fine grid

% X,Y:     1d x & y coordinates which define the whole 2D domain at each level of the V-cycle
% XG, YG:  coordinates of grid points of each level of the V-cycle
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
    z = Z{i};
    xi = XG{i+1};
    yi = YG{i+1};
    zi = ZG{i+1};
    band = BAND{i};
    T = interp3_matrix( x, y, z, xi, yi, zi, p_f2c, band);
   
    TMf2c{i} = T;
end


for i = 1:1:n_level-1
    x = X{i+1};
    y = Y{i+1};
    z = Z{i+1};
    xi = XG{i};
    yi = YG{i};
    zi = ZG{i};
    band = BAND{i+1};

    T = interp3_matrix( x, y, z, xi, yi, zi, p_c2f, band);
    TMc2f{i} = T;

end

end
