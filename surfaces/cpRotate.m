function [cpx, cpy, cpz, varargout] = cpRotate(x,y,z, cpfun, R)
%CPROTATE  Rotate a given closest point function
%   cpf = @(x,y,z) cpCylinder(x,y,z, [-1.2 1])
%   [...] = cpRotate(x,y,z, cpf)
%   [...] = cpRotate(x,y,z, cpf, [th1 th2 th3])
%
%   The three angles represent:
%        th1: rotation CCW around z-axis: Rz
%        th2: rotation CCW around y-axis: Ry
%        th3: rotation CCW around x-axis: Rx
%   Internally, the code builds three rotation matrices Rz Ry and
%   Rx.  The rotation is then performed using R = Rx*Ry*Rz.
%
%   Or you can provide a 3x3 rotation matrix yourself.  Careful: we
%   don't check its a rotation matrix and the closest points will be
%   incorrect if its not!
%   [...] = cpRotate(x,y,z, cpf, R)

  if (nargin == 3)
    cpfun = @cpCylinder;
  end
  if (nargin <= 4)
    angles = [pi/3 pi/10 pi/12];
    R = rotation_matrix(angles);
  elseif (isvector(R))
    R = rotation_matrix(R);
  end

  % apply the inverse (CW) rotation to the points x, then later rotate
  % the closest points by R (CCW).  We right-mult by R which is same
  % as left mult by R'.
  sh = size(x);
  XR = [x(:) y(:) z(:)] * R;
  [cpx, cpy, cpz, varargout{1:nargout-3}] = cpfun(XR(:,1), XR(:,2), XR(:,3));
  CP = [cpx cpy cpz];
  CPR = CP * R';
  cpx = CPR(:,1);
  cpy = CPR(:,2);
  cpz = CPR(:,3);

  cpx = reshape(cpx, sh);
  cpy = reshape(cpy, sh);
  cpz = reshape(cpz, sh);
  for i=1:(nargout-3)
    varargout{i} = reshape(varargout{i}, sh);
  end
