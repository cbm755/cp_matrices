function [cpx, cpy, dist, bdy] = cpSemicircle(x, y, R, xc, yc)
%CPSEMICIRCLE  Closest point function for a semicircle.
%   The semicircle consists of those points with y >= 0.
%   [cpx, cpy, dist, bdy] = cpSemicircle(x, y, R, xc, yc)
%      "bdy" is non-zero for points on the boundary.
%      Inputs R, xc, yc can be omitted and default to (1,0,0)
%
%   Code is vectorized: any size/shape for x should work.


  % radius defaults to 1
  if (nargin < 3)
    R = 1;
  end
  if (nargin == 4)
    error('must specify both xc,yc');
  end
  % center shift, defaults to 0
  if (nargin < 5)
    xc = 0;
    yc = 0;
  end

  % shift to the origin
  x = x - xc;
  y = y - yc;

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

  % note: this way doesn't identify the left and right bdries
  %bdy = (Reg1 | Reg2);
  bdy = Reg1 + 2*Reg2;

  % shift back
  cpx = cpx + xc;
  cpy = cpy + yc;
