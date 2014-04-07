function [cpx, cpy, varargout] = cpRotate2d(x,y, cpfun, R)
%CPROTATE2D  Rotate a given 2D closest point function
%   cpf = @(x,y) cpEllipse(x,y, 1.5, 0.5)
%   [...] = cpRotate2d(x,y, cpf)
%   [...] = cpRotate2d(x,y, cpf, th1)
%   You can specify the angle th1 measured CCW in radians.
%
%   Or you can provide a 2x2 rotation matrix yourself.  Careful: we
%   don't check its a rotation matrix and the closest points will be
%   incorrect if its not!
%
%   [...] = cpRotate2d(x,y, cpf, R)

  if (nargin == 2)
    cpfun = @cpEllipse;
  end
  if (nargin <= 3)
    R = rotation_matrix(pi/3);
  elseif (isscalar(R))
    R = rotation_matrix(R);
  end

  % apply the inverse (CW) rotation to the points x, then later rotate
  % the closest points by R (CCW).  We right-mult by R which is same
  % as left mult by R'.
  sh = size(x);
  XR = [x(:) y(:)] * R;
  [cpx, cpy, varargout{1:nargout-2}] = cpfun(XR(:,1), XR(:,2));
  CP = [cpx cpy];
  CPR = CP * R';
  cpx = CPR(:,1);
  cpy = CPR(:,2);

  cpx = reshape(cpx, sh);
  cpy = reshape(cpy, sh);
  for i=1:(nargout-2)
    varargout{i} = reshape(varargout{i}, sh);
  end
