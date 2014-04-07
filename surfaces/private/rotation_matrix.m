function R = rotation_matrix(angles)
%ROTATION_MATRIX  Return a rotation matrix
%   R = rotation_matrix(pi/3)
%       A rotation in 2D (2x2 matrix) of pi/3 CCW.
%   R = rotation_matrix([pi/3 pi/4 pi/5])
%       A matrix in 3D, here the angles represent:
%          i) rotation CCW around z-axis: Rz
%         ii) rotation CCW around y-axis: Ry
%        iii) rotation CCW around x-axis: Rx
%       R is the R = Rx*Ry*Rz.

  if (isscalar(angles))
    th = angles; c = cos(th); s = sin(th);
    R = [c  -s;  s  c];

  elseif (length(angles == 3))
    th = angles(1); c = cos(th); s = sin(th);
    Rz = [c -s  0;  s  c  0;  0  0  1];
    th = angles(2); c = cos(th); s = sin(th);
    Ry = [c  0 -s;  0  1  0;  s  0  c];
    th = angles(3); c = cos(th); s = sin(th);
    Rx = [1  0  0;  0  c -s;  0  s  c];
    R = Rx*Ry*Rz;

  else
    error('not implemented in dimension higher than 3')
  end

