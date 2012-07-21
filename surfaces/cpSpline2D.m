function [cpx, cpy, dist, varargout] = cpSpline2D(x, y, sp, isclosed, DEBUG)
%CPSPLINE2D  Closest Point function for closed spline in 2D
%   Helper function, see cpEggCurve, cpBeanCurve, cpCosCurve for
%   examples.  Works for both closed and open curves.
%
%   Depends on the spline toolbox
%   Internally, uses cpParamCurve_2D with Newton solves to find cp
%
%   TODO: the open case still uses loops

  if (nargin <= 4)
    DEBUG = 0;
  end

  % first and second derivatives
  sp1 = fnder(sp);
  sp2 = fnder(sp,2);
  %splx = cscvn(pts(:,1)');
  %sply = cscvn(pts(:,1)');

  % TODO: functions with handle references are easier to
  % understand, work for different input sizes and are even
  % slightly faster, so these function handles are deprecated
  S1.type='()';
  S1.subs = {1,':'};
  S2.type='()';
  S2.subs = {2,':'};

  % parameterised curve:
  xs = @(t) subsref(ppval(sp,t), S1);
  ys = @(t) subsref(ppval(sp,t), S2);

  % derivative of parametrisation:
  xp = @(t) subsref(ppval(sp1,t), S1);
  yp = @(t) subsref(ppval(sp1,t), S2);

  % second derivative:
  xpp = @(t) subsref(ppval(sp2,t), S1);
  ypp = @(t) subsref(ppval(sp2,t), S2);

  % start and endpts of parameter variable
  endpt1 = sp.breaks(1);
  endpt2 = sp.breaks(end);

  xv = x(:); yv = y(:);

  if (isclosed)
    useVectorCode = 1;
  else
    useVectorCode = 0;
  end


  if (useVectorCode)
    if (isclosed)
      [cpx, cpy, dist, fail] = ...
          cpParamCurveClosed(xv,yv, ...
              @myxs,@myys,@myxp,@myyp,@myxpp,@myypp, ...
              [endpt1 endpt2], DEBUG);
    else
      error('TODO!')
    end
  else
    % loop and find the CPs
    nx = length(xv);     % number of points

    cpx = zeros(nx,1);
    cpy = zeros(nx,1);
    dist = zeros(nx,1);
    bdy = zeros(nx,1);

    for pt = 1:nx
      if (isclosed)
        [cpx(pt), cpy(pt), dist(pt)] = ...
            cpParamCurveClosed_oldloop(xv(pt),yv(pt), ...
                xs,ys,xp,yp,xpp,ypp, ...
                [endpt1 endpt2], DEBUG);
      else
        [cpx(pt), cpy(pt), dist(pt), bdy(pt)] = ...
            cpParamCurveOpen_oldloop(xv(pt),yv(pt), ...
                xs,ys,xp,yp,xpp,ypp, ...
                [endpt1 endpt2], DEBUG);
      end
    end
  end

  cpx = reshape(cpx, size(x));
  cpy = reshape(cpy, size(x));
  dist = reshape(dist, size(x));

  if (~isclosed)
    bdy = reshape(bdy, size(x));
    varargout = {bdy};
  end

  if (DEBUG >= 1)
    figure(1);
    plot(cpx(~bdy), cpy(~bdy), 'bx');
    hold on;
    axis equal;
    plot(x(~~bdy), y(~~bdy), 'rx');
  end

  function r = myxs(t)
    tt = t(:);
    A = ppval(sp,t);
    r = reshape(A(1,:), size(t));
  end
  function r = myys(t)
    tt = t(:);
    A = ppval(sp,t);
    r = reshape(A(2,:), size(t));
  end
  function r = myxp(t)
    tt = t(:);
    A = ppval(sp1,t);
    r = reshape(A(1,:), size(t));
  end
  function r = myyp(t)
    tt = t(:);
    A = ppval(sp1,t);
    r = reshape(A(2,:), size(t));
  end
  function r = myxpp(t)
    tt = t(:);
    A = ppval(sp2,t);
    r = reshape(A(1,:), size(t));
  end
  function r = myypp(t)
    tt = t(:);
    A = ppval(sp2,t);
    r = reshape(A(2,:), size(t));
  end

end


