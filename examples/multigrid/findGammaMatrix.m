function [GAMMA] = findGammaMatrix(alpha, beta, gamma1, x, y, z, xi, yi, zi, p)

%FINDGAMMAMATRIX: use different gamma's for different closest points,
%hoping to get an M-Matrix in 3D. (fingers crossed.)
% the output GAMMA is a sparse diagonal matrix.

  T1 = cputime();
  dx = x(2)-x(1);   Nx = length(x);
  dy = y(2)-y(1);   Ny = length(y);
  dz = z(2)-z(1);   Nz = length(z);
  ddx = [dx  dy  dz];
  ptL = [x(1) y(1) z(1)];

 

  dim = 3;
  N = p+1;
  % we only compute 8 weights very near the closest point
  EXTSTENSZ = 2^dim;

  %tic
  Ei = repmat((1:length(xi))',1,EXTSTENSZ);
  Ej = zeros(size(Ei));
  gammas = zeros(size(Ei));
  % we need both cubic and linear interp weights
  weights3 = zeros(size(Ei));
  weights1 = zeros(size(Ei));

  %tic
  % this used to be a call to buildInterpWeights but now most of
  % that is done here
  [Ibpt, Xgrid] = findGridInterpBasePt_vec({xi yi zi}, p, ptL, ddx);
  xw3 = LagrangeWeights1D_vec(Xgrid{1}, xi, ddx(1), N);
  yw3 = LagrangeWeights1D_vec(Xgrid{2}, yi, ddx(2), N);
  zw3 = LagrangeWeights1D_vec(Xgrid{3}, zi, ddx(3), N);
  
  xw1 = LagrangeWeights1D_vec(Xgrid{1}+dx, xi, ddx(1), 2);
  yw1 = LagrangeWeights1D_vec(Xgrid{2}+dy, yi, ddx(2), 2);
  zw1 = LagrangeWeights1D_vec(Xgrid{3}+dz, zi, ddx(3), 2);
  
  
  %toc

  %tic
  % this is a good order for memory access: ijk just counts up
  for k=1:2
    for i=1:2
      for j=1:2
        ijk = sub2ind([2,2,2], j, i, k);
        weights3(:,ijk) = xw3(:,i+1) .* yw3(:,j+1) .* zw3(:,k+1);
        weights1(:,ijk) = xw1(:,i) .* yw1(:,j) .* zw1(:,k);
        
        weights_from_L = xw1(:,i).*yw1(:,j) + yw1(:,j).*zw1(:,k) + ...
                         zw1(:,k).*xw1(:,i) - 3*weights1(:,ijk);
                     
        numerator = gamma1*weights1(:,ijk) - alpha*weights_from_L;
        
        %numerator = numerator - beta/4*xw1(:,3-i).*yw1(:,3-j).*zw1(:,3-k);
        
        %numerator = max(numerator,0);
        
        gammas(:,ijk) = numerator ./ weights3(:,ijk);
        
%         ind = abs(weights3(:,ijk)) < 1e-14 & abs(numerator) < 1e-14;
%         gammas(ind,ijk) = 0;
        
        Ej(:,ijk) = ijk;
        
      end
    end
  end
  %toc
  T1 = cputime() - T1;
  %fprintf('done new Ei,Ej,weights, total time: %g\n', T1);

  GAMMA = sparse(Ei(:), Ej(:), gammas(:), length(xi), 8);

  %gamma = min(max(GAMMA,[],2), 5.1) ;
  gamma = max(GAMMA,[],2);
  
  GAMMA = spdiags(gamma/dx^2,0,length(xi),length(xi));
  
    

end
