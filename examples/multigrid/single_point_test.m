function [] = single_point_test(alpha,beta)

    cpx = 0.;
    cpy = 0.;
    cpz = 0.01;
    
    dx = 1;
    x1d = (-2:dx:3)';
    y1d = x1d;
    z1d = x1d;

    nx = length(x1d);
    ny = length(y1d);
    nz = length(z1d);

    [xx,yy,zz] = meshgrid(x1d,y1d,z1d);
    xx = xx(:); yy = yy(:); zz = zz(:);
    band = find ( (xx>=-1) & (xx<=2) & (yy>=-1) & (yy<=2) & (zz>=-1) & (zz<=2) );
    
    xg = xx(band);
    yg = yy(band);
    zg = zz(band);
    
    E1 = interp3_matrix(x1d,y1d,z1d,cpx,cpy,cpz,1,band);
    E3 = interp3_matrix(x1d,y1d,z1d,cpx,cpy,cpz,3,band);
    %L  = laplacian_3d_matrix(x1d,y1d,z1d,2,band);
    L  = laplacian_wider_stencil_3d_matrix(x1d,y1d,z1d,2,alpha,beta,1-alpha-beta,band,band);
    
    for lambda = 2:0.1:6
        tmp = E1*L + lambda*E3;
        ind = tmp<-1e-10;
        [lambda, tmp(ind)]
        [xg(ind), yg(ind), zg(ind)]
        %tmp(ind)
        if (nnz(tmp<-1e-10) == 0)
            lambda
            break
        end
    end
    
end