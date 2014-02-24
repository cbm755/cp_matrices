function [GAMMA1, GAMMA2] = findGammaMatrix2(x, y, z, xi, yi, zi, p)

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
                     
        numerator = 6*weights1(:,ijk); %- weights_from_L;
        
        %numerator = max(numerator,0);
        
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
  hx = xi - Xgrid{1} - dx;
  hy = yi - Xgrid{2} - dx;
  hz = zi - Xgrid{3} - dx;
  
  % weights from diagonal entries of Laplacian.
  diagwL = 2 * ( 1./((dx+hx).*(2*dx-hx)) + 1./((dy+hy).*(2*dy-hy)) + ...
                1./((dz+hz).*(2*dz-hz)) );
            
  % weights from other 6 directions.
  xwL = cell(2,1); ywL = cell(2,1); zwL = cell(2,1);
  xwL{1} = 2/3/dx./(dx+hx);            % negative x-direction
  xwL{2} = 2/3/dx./(2*dx-hx);          % postive x-direction
  ywL{1} = 2/3/dy./(dy+hy);            % negative y-direction
  ywL{2} = 2/3/dy./(2*dy-hy);          % postive y-direction
  zwL{1} = 2/3/dz./(dz+hz);            % negative z-direction
  zwL{2} = 2/3/dz./(2*dz-hz);          % postive z-direction

  % compute the actual weights of the Laplacian at CP along the x-direction
  for cnt = 1:2
      % we want i = 1 & 4 correspondingly
      i = 3*cnt - 2;
      for k = [1,4]
          for j = [1,4]
              ijk = sub2ind([N,N,N], j, i, k);
              weightsL(:,ijk) = weightsL(:,ijk) + xwL{cnt}(:) .* yw3(:,j) .* zw3(:,k);
          end
      end
  end
  
  % compute the actual weights of the Laplacian at CP along the y-direction
  for cnt = 1:2
      % we want j = 1 & 4 correspondingly
      j = 3*cnt - 2;
      for k = [1,4]
          for i = [1,4]
              ijk = sub2ind([N,N,N], j, i, k);
              weightsL(:,ijk) = weightsL(:,ijk) + ywL{cnt}(:) .* xw3(:,i) .* zw3(:,k);
          end
      end
  end
  
  % compute the actual weights of the Laplacian at CP along the z-direction
  for cnt = 1:2
      % we want k = 1 & 4 correspondingly
      k = 3*cnt - 2;
      for i = 1:N
          for j = 1:N
              ijk = sub2ind([N,N,N], j, i, k);
              weightsL(:,ijk) = weightsL(:,ijk) + zwL{cnt}(:) .* xw3(:,i) .* yw3(:,j);
          end
      end
  end
  %toc
  
  % extract the 8 corner weights 
  weightsL_corner = weightsL(:,[1,4,13,16,49,52,61,64]);
  
  DiagwLinverse = spdiags([1./diagwL],0,length(xi),length(xi));
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
