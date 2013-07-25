function [cpx, cpy, d, bdy, ls] = cpYinYang(x, y, R, cen)
%CPYINYANG  Closest Point function for a mixed-codim Yin Yang.
%   [cpx, cpy, dist, bdy, levelset] = cpYinYang(x, y)
%      A unit yin yang centered at the origin
%   [...] = cpCircle(x, y, R)
%   [...] = cpCircle(x, y, R, Cen)
%      A yin yang with radius R centered at Cen = [xc,yc].
%
%   Outputs:
%     `dist' is unsigned distance.  Note zero in the "solid" part.
%     `bdy' labels different components.
%     `levelset' is a signed distance representation.  Note its
%        nonzero (and negative) inside the "solid" part and thus
%        different than dist.
%

  % defaults
  if (nargin < 3)
    R = 1;
  end
  if (nargin < 4)
    cen = [0,0];
  end

  % shift to the origins
  x = x - cen(1);
  y = y - cen(2);

  [cpx1,cpy1,sd1] = cpCircle(x, y, R);
  [cpx2,cpy2,sd2] = cpCircle(x, y, R/2, [R/2 0]);
  [cpx3,cpy3,sd3] = cpCircle(x, y, R/2, [-R/2 0]);

  cpx = zeros(size(x));
  cpy = zeros(size(x));
  % for debugging
  d = -1000*ones(size(x));
  ls = -1000*ones(size(x));
  bdy = -1000*ones(size(x));

  % case 1: outside big circle 1, cp is cpx1
  I = (sd1 >= 0);
  cpx(I) = cpx1(I);
  cpy(I) = cpy1(I);
  d(I) = sd1(I);
  ls(I) = sd1(I);
  % two subcases to get different labels
  bdy(I) = 1;
  bdy(I & (y <= 0)) = 3;

  % next two are inside the black yin
  I = (sd1 < 0) & (sd3 <= 0) & (y < 0);
  cpx(I) = x(I);
  cpy(I) = y(I);
  d(I) = 0;
  ls(I) = sd3(I);
  bdy(I) = 0;

  I = (sd1 < 0) & (sd2 > 0) & (y >= 0);
  cpx(I) = x(I);
  cpy(I) = y(I);
  d(I) = 0;
  bdy(I) = 0;
  % level set is setup below for this case


  %% next three are the white yang
  % inside the yang and in top half-plane
  I = (sd2 <= 0) & (y >= 0);
  cpx(I) = cpx2(I);
  cpy(I) = cpy2(I);
  d(I) = -sd2(I);
  ls(I) = -sd2(I);
  bdy(I) = 2;

  % inside the yang, bottom half-plane and closest to the yin
  I = (sd1 < 0) & (y < 0) & (sd3 > 0) & (sd3 <= -sd1);
  cpx(I) = cpx3(I);
  cpy(I) = cpy3(I);
  d(I) = sd3(I);
  ls(I) = sd3(I);
  bdy(I) = 2;

  % inside the yang, bottom half-plane and closest to the big circle
  I = (sd1 < 0) & (y < 0) & (sd3 > 0) & (sd3 > -sd1);
  cpx(I) = cpx1(I);
  cpy(I) = cpy1(I);
  d(I) = -sd1(I);
  ls(I) = -sd1(I);
  bdy(I) = 3;

  % few extra cases to get signed distance function
  I = (sd1 < 0) & (y >= 0) & (sd2 > 0) & (sd2 <= -sd1);
  ls(I) = -sd2(I);

  I = (sd1 < 0) & (y >= 0) & (sd2 > 0) & (sd2 > -sd1);
  ls(I) = sd1(I);

  % sanity checks in case we missed some cases
  if (any(any(d < 0)))
    error('negative distance');
  end
  if (any(any(ls < -999)))
    error('level set too small');
  end
  if (any(any(bdy < 0)))
    error('invalid bdy labelling');
  end

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
