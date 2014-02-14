function [TMf2c, TMc2f] = helper_set_TM3d(X, Y, Z, XCP, YCP, ZCP, BAND, BDYG, p_f2c, p_c2f)

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
    
    %T = interp3_matrix_band( x, y, z, xi, yi, zi, p, band, false );
    T = interp3_matrix( X{i}, Y{i}, Z{i}, XCP{i+1}, YCP{i+1}, ZCP{i+1}, p_f2c);
    %T(BDYG{i+1},:) = -T(BDYG{i+1},:);

    TMf2c{i} = T(:,BAND{i});
    clear T
end


for i = 1:1:n_level-1

    %T = interp3_matrix_band( x, y, z, xi, yi, zi, p_mg, band, false );
    T = interp3_matrix( X{i+1}, Y{i+1}, Z{i+1}, XCP{i}, YCP{i}, ZCP{i}, p_c2f);
    %T(BDYG{i},:) = -T(BDYG{i},:);
    TMc2f{i} = T(:,BAND{i+1});
    clear T
end

end
