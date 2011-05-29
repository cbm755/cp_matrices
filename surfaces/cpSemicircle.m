function [cpx, cpy, dist, bdy] = cpSemicircle(x, y, R)
%CPSEMICIRCLE  Closest point function for a semicircle.
%   The semicircle consists of those points with y >= 0.
%   "bdy" is non-zero for points on the boundary.

  [cpx, cpy] = cpCircle(x, y, R);

  % should not include the endpoints (hence y<0, not y<=0), this is b/c
  % cpbar would then be a boundary point again.  But probably doesn't
  % matter much in practice.
  Reg1 = (x < 0) & (y < 0);
  Reg2 = (x >= 0) & (y < 0);

  cpx(Reg1) = -R;
  cpy(Reg1) = 0;

  cpx(Reg2) = R;
  cpy(Reg2) = 0;

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 );

  bdy = (Reg1 | Reg2);
