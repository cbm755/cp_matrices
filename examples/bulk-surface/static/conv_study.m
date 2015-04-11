function [err_inf, err_L2, ratio_inf, ratio_L2, ddx] = conv_study(f, n_level, dx_coarsest)

    % e.g.: [err, ratio, ddx] = conv_study(@poisson_disk_eg1, 10);
    % [err, ratio, ddx] = conv_study(@poisson_disk_eg2, 10, 0.05);
    
    err_inf = [];
    err_L2 = [];
    ddx = [];
    if nargin < 3
        dx_coarsest = 0.1;
    end
    
    for i = 1:1:n_level
        dx = dx_coarsest/2^(i-1);
        [tmp_inf, tmp_L2] = f(dx);
        err_inf = [err_inf,tmp_inf];
        err_L2 = [err_L2,tmp_L2];
        ddx = [ddx,dx];
    end
    
    ratio_inf = err_inf(1:end-1) ./ err_inf(2:end);
    ratio_L2 = err_L2(1:end-1) ./ err_L2(2:end);
    
end
