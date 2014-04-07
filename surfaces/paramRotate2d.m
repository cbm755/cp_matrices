function [x, y] = paramRotate2d(N, paramfun, R)
%PARAMROTATE2D  Rotate a given parametered curve
%   paramf = @(n) paramEllipse(n, 1.5, 0.5)
%   [x,y] = paramRotate2d(n, paramf)

  if (nargin == 1)
    paramfun = @paramEllipse;
  end
  if (nargin <= 2)
    R = rotation_matrix(pi/3);
  elseif (isscalar(R))
    R = rotation_matrix(R);
  end

  [x, y] = paramfun(N);
  sh = size(x);

  % we want to left-mult R on column vector x.  But easiest to put
  % each point as a row of X.   (R*x)' = x' * R' = X * R.
  X = [x(:) y(:)];
  XR = X * R';
  x = XR(:,1);  y = XR(:,2);
  x = reshape(x, sh);
  y = reshape(y, sh);

