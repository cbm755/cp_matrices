function [cpx, cpy, sdist, tt] = cpEllipse_nonEuclidean(x, y, a, b, cen)
%CPELLIPSE_NONEUCLIDEAN Closest points for ellipse, w/ non-Euclidean distance
%   [cpx, cpy, sdist] = cpEllipse_nonEuclidean(x, y, a, b)
%      An ellipse centered at the origin with major axis 'a' and
%      minor axis 'b'.  'a' and 'b' default to 1.5 and 0.75.
%   [cpx, cpy, sdist] = cpEllipse_nonEuclidean(x, y, a, b, cen)
%      An ellipse centered at 'cen'.
%
%   Note: returns signed non-Euclidean distance (with negative inside).

  % defaults
  if (nargin < 4)
    if (nargin == 3)
      error('must provide both or either of a,b');
    end
    a = 1.5;
    b = 0.75;
  end
  if (nargin < 5)
    cen = [0,0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);

  sdist = sqrt(x.^2/a^2 + y.^2/b^2) - 1;

  % this function is zero when t is the parameter of the closest point
  g = @(x,y,t) (a*cos(t)).^(a^2) .* y.^(b^2)  ...
             - (b*sin(t)).^(b^2) .* x.^(a^2);

  %% first, put everything in the first quadrant
  % todo: on axis?
  Q1 = (x >= 0) & (y >= 0);
  Q2 = (x < 0)  & (y >= 0);
  Q3 = (x < 0)  & (y < 0);
  Q4 = (x >= 0) & (y < 0);
  Iax = (y == 0);
  Iay = (x == 0);

  x(Q2) = -x(Q2);
  y(Q4) = -y(Q4);
  x(Q3) = -x(Q3);
  y(Q3) = -y(Q3);

  assert(all(all(x >= 0)))
  assert(all(all(y >= 0)))

  % deal with axis exactly
  cpx = zeros(size(x));
  cpy = zeros(size(x));
  tt = zeros(size(x));
  cpx(Iax) = a;
  cpy(Iax) = 0;
  cpx(Iay) = 0;
  cpy(Iay) = b;
  tt(Iax) = 0;
  tt(Iay) = pi/2;

  % get the ones that aren't on axis, we will do bisection on these
  X = x( ~Iax & ~Iay );
  Y = y( ~Iax & ~Iay );

  if ~isempty(X)
    time_bisect = cputime();
    t1 = zeros(size(X));
    t2 = pi/2*ones(size(X));
    g1 = g(X,Y,t1);
    g2 = g(X,Y,t2);

    % if this fails, may need to consider hacks for end points
    assert(all(g1 .* g2 < 0));

    n = 0;
    while (1)
      if max(abs(t1-t2)) < 2*eps
        fprintf('Bisection: converged in n = %d iterations, %g sec\n', ...
                n, cputime()-time_bisect)
        t = (t1+t2)/2;
        break
      end
      tnew = (t1 + t2)/2;
      gnew = g(X,Y,tnew);
      I0 = gnew == 0;
      I1 = g1 .* gnew < 0;
      I2 = g2 .* gnew < 0;
      % close the interval on exact values
      g1(I0) = gnew(I0);
      g2(I0) = gnew(I0);
      t1(I0) = tnew(I0);
      t2(I0) = tnew(I0);
      % usual bisection
      g1(I2) = gnew(I2);
      t1(I2) = tnew(I2);
      g2(I1) = gnew(I1);
      t2(I1) = tnew(I1);
      n = n + 1;
    end

    %max(abs(g(X,Y,t)))
    assert(max(abs(g(X,Y,t))) < 1e4*eps)

    cpx( ~Iax & ~Iay ) = a*cos(t);
    cpy( ~Iax & ~Iay ) = b*sin(t);
    tt( ~Iax & ~Iay ) = t;
  end

  cpx(Q2) = -cpx(Q2);
  cpy(Q4) = -cpy(Q4);
  cpx(Q3) = -cpx(Q3);
  cpy(Q3) = -cpy(Q3);

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);

