function [cpx, cpy, sdist] = cpEllipse(x, y, a, b, cen)
%CPELLIPSE  Closest Point function for an ellipse.
%   [cpx, cpy, sdist] = cpEllipse(x, y, a, b)
%      An ellipse centered at the origin with major axis 'a' and
%      minor axis 'b'.  'a' and 'b' default to 1.5 and 0.75.
%   [cpx, cpy, sdist] = cpEllipse(x, y, a, b, cen)
%      An ellipse centered at 'cen'.
%
%   Internally, uses cpParamCurve with Newton solves to find cp
%
%   Note: returns signed distance (with negative inside).

  % defaults
  if (nargin < 4)
    if (nargin == 3)
      error('must provide both a and b (or neither)');
    end
    a = 1.5;
    b = 0.75;
  end
  if (nargin < 5)
    cen = [0 0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);

  % parameterised curve:
  xs = @(t) a*cos(t);
  ys = @(t) b*sin(t);

  % derivative of parametrisation:
  xp = @(t) -a*sin(t);
  yp = @(t) b*cos(t);

  % second derivative:
  xpp = @(t) -a*cos(t);
  ypp = @(t) -b*sin(t);

  % start and endpts of parameter variable
  endpt1 = 0;
  endpt2 = 2*pi;

  xv = x(:); yv = y(:);


  useVectorCode = 1;

  if (useVectorCode)
    [cpx, cpy, dist, fail] = ...
        cpParamCurveClosed(xv,yv, xs,ys,xp,yp,xpp,ypp, [endpt1 endpt2]);

  else
    warning('this code is very slow');

    nx = length(xv);     % number of points
    for pt = 1:nx
      [cpx(pt), cpy(pt), dist(pt), fail(pt)] = ...
          cpParamCurveClosed_oldloop(xv(pt),yv(pt),xs,ys,xp,yp,xpp,ypp, [endpt1 endpt2]);
      %[t1,t2,t3,t4] = cpParamCurve_2D(x1d(pt),y1d(pt),xs,ys,xp,yp,xpp,ypp,endpt1,endpt2,1);
      %cpx(pt) = t1;
      %cpy(pt) = t2;
      %dist(pt) = t3;
      %fail(pt) = t4;
    end
  end

  %max(abs(cpx-cpx2))
  %max(abs(cpy-cpy2))
  %max(abs(dist-dist2))
  %keyboard

  cpx = reshape(cpx, size(x));
  cpy = reshape(cpy, size(x));
  dist = reshape(dist, size(x));
  fail = reshape(fail, size(x));

  % in the top (above 45degree lines), y vs cpy will be a good indicator
  % of sign
  sdist = ...
      (y>=x & y>-x).*(   sign((y-cpy)).*dist  ) + ...   % top
      (y>x  & y<=-x).*( -sign((x-cpx)).*dist  ) + ...   % right
      (y<=x & y<-x).*(  -sign((y-cpy)).*dist  ) + ...   % bottom
      (y<x  & y>=-x).*(  sign((x-cpx)).*dist  );        % left

  % above cases miss the origin
  wh = (x==0 & y==0);
  %wh = (x.^2 + y.^2 <= (1/8*min(a,b))^2);   %use small circle instead
  sdist = (~wh) .* sdist + (wh) .* (-dist);


  if (any(any(fail)))
    warning('some cpEllipse points failed');
    figure(1); clf;
    plot(cpx(~fail), cpy(~fail), 'bx');
    hold on;
    axis equal;
    plot(x(~~fail), y(~~fail), 'rx');
    %keyboard
  end

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
