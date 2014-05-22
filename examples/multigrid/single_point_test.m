function [] = single_point_test(alpha,beta)

    cpx = 0.1;
    cpy = 0.9;
    cpz = 0.2;
    
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
    
    [Eaverage1,DiagwLinverse1] = buildCPmatrixFromLaplacian3d2ndOrder3rdAttempt(x1d,y1d,z1d,cpx,cpy,cpz,3,band);
    [GAMMA1, GAMMA2] = find2GammaMatrices(x1d, y1d, z1d, cpx, cpy, cpz, 3);
    
%    [Eaverage1,DiagwLinverse1] = buildCPmatrixFromLaplacian3d2ndOrder(x1d,y1d,z1d,cpx,cpy,cpz,3,band);
%    [GAMMA1, GAMMA2] = find2GammaMatrices2ndOrder(x1d, y1d, z1d, cpx, cpy, cpz, 3);

%    [Eaverage2,DiagwLinverse2] = buildCPmatrixFromLaplacian3d2ndOrder2ndAttempt(x1d,y1d,z1d,cpx,cpy,cpz,3,band);
%    [GAMMA1, GAMMA2] = find2GammaMatrices2ndOrder2ndAttempt(x1d, y1d, z1d, cpx, cpy, cpz, 3);
%    Lcp = buildLaplacianMatrixAtCP(x1d,y1d,z1d,cpx,cpy,cpz,3,band);


    GAMMA1
    GAMMA2
    
    for lambda1 = 6:0.1:6
%      lambda1
      for lambda2 = 1:0.1:1
%        for lambda3 = 0:0.1:1
          %tmp = E1*L + lambda1*E3 + lambda2*(DiagwLinverse1*Eaverage1) - lambda3*DiagwLinverse1*Lcp;
          tmp = E1*L + GAMMA1*E3 +GAMMA2*DiagwLinverse1*Eaverage1;
          reshape(full(tmp),4,4,4)
          %tmp = (1-lambda2)*E1*L + lambda2*Lcp + lambda1*E3;
          ind = tmp<-1e-14;
%           nnz(ind)
%           [lambda1, w]
%           min(tmp(ind))
           [xg(ind), yg(ind), zg(ind)]
%           tmp(ind)
          if (nnz(tmp<-1e-14) == 0)
            
            %[lambda1,lambda2]
            disp('good')
          end
%        end
      end
    end
    
end