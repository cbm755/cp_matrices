function [cpx, cpy, d, bdy, ls] = cpMixedCodim2d(x, y, R, a, b, cen)
%CPMIXEDCODIM2D  A mixed codimension test case

  % defaults
  if (nargin == 2)
    R = 1;
    a = 0.9;
    b = 0.6;
  elseif (nargin > 2) && (nargin < 5)
    error('must provide none or all of R, a, b');
  end
  if (nargin < 6)
    cen = [0 0];
  end

  % shift to the origins
  x = x - cen(1);
  y = y - cen(2);

  % 3 conditions: agreement of the two curves equations and then that
  % the tangents be parallel (normal of circle dotted with tangent of
  % ellipse).  Put this into a CAS.
  A = (R^2-a^2) / (a^2 - b^2);
  B = (a^4 - b^2*R^2) / (a^2 - b^2);

  if (A < 0) || (B < 0)
    error('invalid R,a,b');
  end
  L2 = sqrt(A);  % +/- root controls up/down
  L3 = sqrt(B);  % +/- root controls left/right

  % 0 = b^2*r^2-a^4 + (-b^2+a^2)*Z^2
  theta = atan2(L2*a/R,L3/R);
  h = L2*(-b^2+a^2)/a;
  %h3 = sqrt( (R^2-a^2)*(a^2-b^2) ) / a
  beta = atan2(L2/a*b, L3/a);

  [cpx1,cpy1,sd1] = cpCircle(x, y, R);
  [cpx2,cpy2,sd2] = cpEllipse(x, y, a, b, [0 h]);

  yT = R*sin(theta);  % threshold

  cpx = zeros(size(x));
  cpy = zeros(size(x));
  % for debugging
  d = -1000*ones(size(x));
  ls = -1000*ones(size(x));
  bdy = -1000*ones(size(x));

  %% divide into several cases
  % outside big circle 1 and below the threshold y, cp is cpx1
  I = (sd1 >= 0) & (cpy1 < yT);
  cpx(I) = cpx1(I);
  cpy(I) = cpy1(I);
  d(I) = sd1(I);
  ls(I) = sd1(I);
  bdy(I) = 2;

  % todo: gets some extra points, to be overwritten later??
  %I = (sd1 >= 0) & (cpy1 >= yT);
  I = (cpy2 >= yT);
  cpx(I) = cpx2(I);
  cpy(I) = cpy2(I);
  d(I) = sd2(I);
  ls(I) = sd2(I);
  bdy(I) = 1;

  % case 2: inside circle outside ellipse, below threshold
  I = (sd1 < 0) & (sd2 >= 0) & (cpy2 < yT);
  cpx(I) = x(I);
  cpy(I) = y(I);
  d(I) = 0;
  ls(I) = -min(-sd1(I), sd2(I));
  bdy(I) = 0;

  % inside ellipse, top part: done except sign on the distance function
  I = (sd2 < 0) & (cpy2 >= yT);
  d(I) = -sd2(I);
  ls(I) = -sd2(I);

  % inside ellipse, bottom part
  I = (sd2 < 0) & (cpy2 < yT);
  cpx(I) = cpx2(I);
  cpy(I) = cpy2(I);
  d(I) = -sd2(I);
  ls(I) = -sd2(I);
  bdy(I) = 3;

  % inside ellipse: done except sign on the level set
  %I = (sd2 < 0);
  %ls(I) = -sd2(I);


  % sanity checks in case we missed some cases
  if (any(any(d < 0)))
    warning('negative distance');
  end
  if (any(any(ls < -999)))
    warning('level set too small');
  end
  if (any(any(bdy < 0)))
    warning('invalid bdy labelling');
  end



  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
