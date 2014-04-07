function [x, y, z] = paramRotate(N, paramfun, R)
%PARAMROTATE  Rotate a given parametered surface
%   paramf = @(n) paramCylinder(n, [-1.2 1])
%   [x,y,z] = paramRotate(n, paramf)


  if (nargin == 1)
    paramfun = @paramCylinder;
  end
  if (nargin <= 2)
    angles = [pi/3 pi/10 pi/12];
    R = rotation_matrix(angles);
  elseif (isvector(R))
    R = rotation_matrix(R);
  end

  [x, y, z] = paramfun(N);
  sh = size(x);

  % we want to left-mult R on column vector x.  But easiest to put
  % each point as a row of X.   (R*x)' = x' * R' = X * R.
  X = [x(:) y(:) z(:)];
  XR = X * R';
  x = XR(:,1);  y = XR(:,2);  z = XR(:,3);
  x = reshape(x, sh);
  y = reshape(y, sh);
  z = reshape(z, sh);

