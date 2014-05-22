function [GAMMA1, GAMMA2] = find2GammaMatrices2ndOrder2ndAttempt(x, y, z, xi, yi, zi, p)

%FINDGAMMAMATRIX: use different gamma's for different closest points,
%need two gamma matrices for the new attempt.
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
  gammas1 = zeros(size(Ei));
  gammas2 = zeros(size(Ei));
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

  %% Compute GAMMA1 (s.t. the center 8 points have positive coeffs)
  %tic
  % this is a good order for memory access: ijk just counts up
  for k=1:2
    for i=1:2
      for j=1:2
        ijk = sub2ind([2,2,2], j, i, k);
        weights3(:,ijk) = xw3(:,i+1) .* yw3(:,j+1) .* zw3(:,k+1);
        weights1(:,ijk) = xw1(:,i) .* yw1(:,j) .* zw1(:,k);
        
%         weights_from_L = xw1(:,i).*yw1(:,j) + yw1(:,j).*zw1(:,k) + ...
%                          zw1(:,k).*xw1(:,i) - 3*weights1(:,ijk);
                     
        numerator = 6*weights1(:,ijk);% - weights_from_L;
        
%         numerator = max(numerator,0);
        
        gammas1(:,ijk) = numerator ./ weights3(:,ijk);
        
%         ind = abs(weights3(:,ijk)) < 1e-14 & abs(numerator) < 1e-14;
%         gammas(ind,ijk) = 0;
        
        Ej(:,ijk) = ijk;
        
      end
    end
  end
  %toc
  T1 = cputime() - T1;
  %fprintf('done new Ei,Ej,weights, total time: %g\n', T1);

  GAMMA1 = sparse(Ei(:), Ej(:), gammas1(:), length(xi), 8);

  %gamma = min(max(GAMMA,[],2), 5.1) ;
  gamma1 = max(GAMMA1,[],2);
  %gamma1(:) = 6;
  
  GAMMA1 = spdiags(gamma1/dx^2,0,length(xi),length(xi));
  
  %% Compute GAMMA2 (s.t. the corner 8 points have positive coeffs)
  weightsL = zeros(length(xi),64);
  % relative coordinates w.r.t. the central block.
  cpxn = cell(6,1); 
  cpyn = cell(6,1);
  cpzn = cell(6,1);
  for i = 1:6
      cpxn{i} = xi;
      cpyn{i} = yi;
      cpzn{i} = zi;
  end
  cpxn{1} = xi - dx;
  cpxn{2} = xi + dx;
  cpyn{3} = yi - dy;
  cpyn{4} = yi + dy;
  cpzn{5} = zi - dz;
  cpzn{6} = zi + dz;
    
  %tic
  % compute the weights and positions
  for cnt = 1:6
      xw = LagrangeWeights1D_vec(Xgrid{1}, cpxn{cnt}, ddx(1),N);
      yw = LagrangeWeights1D_vec(Xgrid{2}, cpyn{cnt}, ddx(2),N);
      zw = LagrangeWeights1D_vec(Xgrid{3}, cpzn{cnt}, ddx(3),N);
      for k=1:N
          for i=1:N
              for j=1:N
                  ijk = sub2ind([N,N,N], j, i, k);
                  weightsL(:,ijk) = weightsL(:,ijk) + xw(:,i) .* yw(:,j) .* zw(:,k);
              end
          end
      end
  end
  weightsL = weightsL/dx^2;
  % extract the 8 corner weights 
  weightsL_corner = weightsL(:,[1,4,13,16,49,52,61,64]);
  
  DiagwLinverse = dx^2/6*speye(length(xi),length(xi));
  weightsL_corner = DiagwLinverse*weightsL_corner;
  
  weights3_corner = zeros(size(weightsL_corner));
 
  for cnt_k = 1:2
      k = 3*cnt_k - 2;
      for cnt_i = 1:2
          i = 3*cnt_i - 2;
          for cnt_j = 1:2
              j = 3*cnt_j - 2;
              ijk = sub2ind([2,2,2], cnt_j,cnt_i,cnt_k);
              weights3_corner(:,ijk) = xw3(:,i).*yw3(:,j).*zw3(:,k);
              gammas2(:,ijk) = - weights3_corner(:,ijk) ./ weightsL_corner(:,ijk);
              ind = abs(weightsL_corner(:,ijk))<1e-14;
              gammas2(ind, ijk) = 0;
          end
      end
  end
  
  % notice that GAMMA1 already contains a factor of dx^2
  gammas2 = GAMMA1*gammas2;
  GAMMA2 = sparse(Ei(:), Ej(:), gammas2(:), length(xi), 8);

  gamma2 = max(GAMMA2,[],2);
  
  % no need to divide by dx^2 since GAMMA1 already contains that
  GAMMA2 = spdiags(gamma2,0,length(xi),length(xi));            

end
