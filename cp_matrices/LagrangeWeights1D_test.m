function w = LagrangeWeights1D_test(xg, x, dx, N)
%LAGRANGEWEIGHTS1D   Barycentric Lagrange interpolation weights
%   w = LagrangeWeights1D(xg, x, dx, N)
%   1D Barycentric Lagrange interpolation weights.
%   xg is the base grid point, the left end grid point
%     (see also "B" in the diagram in findGridInterpBasePt.m)
%   x is the point of interpolation
%   dx is the grid spacing (uniform)
%   N is the number of points in the stencil (degree + 1)

  w = zeros(size(x,1),N);
  flag = true(size(x,1),1);
  vec = (0:N-1);
  dif = x(:,ones(N,1)) - (xg(:,ones(N,1)) + vec(ones(size(x)),:)*dx);
  [i j] = find( dif == 0 );
  w(i,j+1) = 1;
  flag(i) = false;
  % Is exactly on a grid point, then return binary weights
%   for j=0:(N-1)
%     if (x == (xg+j*dx))
%       w = zeros(1,N);
%       w(j+1) = 1;
%       return
%     end
%   end

  % I think its faster to hardcode than use nchoosek

  %# for i in range(0,N):
  %  #    w[i] = (-1)**i * comb(N-1,i)
  if     (N == 4), w1 = [1, -3, 3, -1];
  elseif (N == 5), w1 = [1, -4, 6, -4, 1];
  elseif (N == 1), w1 = 1;
  elseif (N == 2), w1 = [1, -1];
  elseif (N == 3), w1 = [1, -2, 1];
  elseif (N == 6), w1 = [1, -5, 10, -10, 5, -1];
  elseif (N == 7), w1 = [1, -6, 15, -20,  15,  -6,   1];
  elseif (N == 8), w1 = [1, -7, 21, -35,  35, -21,   7,  -1];
  else error('need to hardcode more weights');
  end

  %for j=0:(N-1)
  %  w(j+1) = w(j+1) / ( x - (xg + j*dx) );
  %end
  w1 = repmat(w1,size(w,1),1);
  w(flag,:) = w1(flag,:);
  w(flag,:) = w(flag,:) ./ dif(flag,:);
  w(flag,:) = w(flag,:) ./ repmat(sum(w(flag,:),2),1,N);
  return

end

