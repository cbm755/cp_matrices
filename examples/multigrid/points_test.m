function [] = points_test()

%     h = 1e-3;
%     x1d = 0;
%     y1d = 0:h:1;
%     z1d = 0:h:1;
%     [xx, yy, zz] = meshgrid(x1d,y1d,z1d);
%     
%     cpx = xx(:);
%     cpy = yy(:);
%     cpz = zz(:);
    
    n = 1e5;
    cpx = 1e-7*rand(n,1);
    cpy = 1e-4*rand(n,1);
    cpz = 1e-2*rand(n,1);
    
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
    L  = laplacian_3d_matrix(x1d,y1d,z1d,2,band);
    
    Eaverage = buildCPmatrixFromLaplacian3d(x1d,y1d,z1d,cpx,cpy,cpz,3,band);
    
    [GAMMA1, GAMMA2] = findGammaMatrix2(x1d, y1d, z1d, cpx, cpy, cpz, 3);

    tmp = E1*L + GAMMA1*E3 +GAMMA2*Eaverage;
    ind = tmp<-1e-14;
    min(tmp(ind))
    [xg(ind), yg(ind), zg(ind)]
    if (nnz(tmp<-1e-14) == 0)
        disp('good')
    end
      
    
end