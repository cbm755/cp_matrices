function [cpx, cpy, cpz, dist, s, fail] = cpParam3DCurveClosed(xx, yy, zz, paramf, endpt, DEBUG)
%CPPARAM3DCURVECLOSED   Find closest points of parameterised curves
%   Uses Newton's method.

% TODO: document!
%
% TODO: this is based on Ingrid's ParamCurve in 2D.  Easily
% genealizes to n-d.  Maybe should write a single function that
% does the n-d version to replace these...
%
% endpt = [pt1 pt2]: endvalues of parameter (e.g., [0 2*pi])
%

  if (nargin < 6)
    DEBUG = 0;
  end

  if (any(size(xx) ~= [1 1]))
    error('not vectorized');
  end

  xpt = [xx; yy; zz];

  % distance squared from (x,y) to curve
  %d2 = @(t,x,y) (xs(t) - x).^2 + (ys(t) - y).^2;

  % derivative of dist wrt t
  %g = @(t,x,y) 2*(xs(t) - x).*xp(t) + 2*(ys(t) - y).*yp(t);

  % second derivative:
  %gp = @(t,x,y) 2*xp(t).*xp(t) + 2*(xs(t) - x).*xpp(t) + 2*yp(t).*yp(t) + 2*(ys(t) - y).*ypp(t);

  ss = linspace(endpt(1), endpt(2), 500);
  [x,xu,xuu] = paramf(ss);
  %dd = sum( (x-xpt).^2 );
  % vectorized so need this: but this is dimensional dependent
  dd = (x(1,:) - xpt(1)).^2 + ...
       (x(2,:) - xpt(2)).^2 + ...
       (x(3,:) - xpt(3)).^2;
  %dd = d2(ss, x, y);
  [mindd_guess, i] = min(dd);
  s_guess = ss(i);


  %% Newton's method
  tol = 1e-14;
  maxn = 100;

  n = 0;
  s = s_guess;

  while (true)
    n = n + 1;

    [x, xp, xpp] = paramf(s);

    f = 2*(x-xpt)' * xp;
    fp = 2*(x-xpt)' * xpp + 2*xp' * xp;

    if (abs(fp) > tol)
      snew = s - f / fp;
    else
      % second deriv is zero
      snew = s;
    end

    %% Force parameter to be periodic
    % TODO: Newton's method if probably not robust if the curve is not
    % smooth at this point.  Needed if the functions of the
    % parameterization cannot be evaluated outside of endpt (i.e.,
    % if they are not periodoc)
    %snew2 = argPeriodic(snew-endpt(1), endpt(2)-endpt(1)) + endpt(1);
    %[s snew s-snew snew2 snew2-snew]
    %snew = snew2;

    if (DEBUG >= 10)
      figure(10);
      plot([xs(s)], [ys(s)],'bo');
      drawnow();
    end

    if (n > maxn)
      fail = 1;
      fprintf('too many iterations: (n,snew,x,y) = %d, %f, %f, %f \n', n, snew, x, y);
      f
      fp
      warning('max iterations in Newton solve: CP is likely wrong!');
      keyboard
      break;
    elseif (abs(s-snew) < tol)
      fail = 0;
      %fprintf('converged: (cpx,cpy,s,n)= %f, %f, %f, %d \n',xs(snew),ys(snew),snew,n);
      break;
    end
    % update
    s = snew;
  end

  % TODO: if we had a list of parameter values corresponding to
  % non-smooth points (i.e., "breaks" from a spline representation), we
  % could loop over them and check if any are closer.  (Newton's
  % method will likely have trouble).

  %cpx = xs(s); cpy = ys(s);
  cp = paramf(s);
  %cpx = cp(1);  cpy = cp(2);
  dd = sum( (cp - xpt).^2 );
  %dist = sqrt((cpx - xx).^2 + (cpy - yy).^2);

  %if (fail)
  %  varargout = {1};
  %else
  %  varargout = {0};
  %end

  cpx = cp(1); cpy = cp(2); cpz = cp(3);
  dist = sqrt(dd);

  % TODO
  if (DEBUG <= 0)
  assert(mindd_guess+tol >= dd, ...
         'initial guess was better: %g, (x,y,z)=(%g,%g,%g), s=%g\n', ...
         mindd_guess - dd, xx, yy,zz,  s);
  else
  if ~(mindd_guess+tol >= dd)
    fprintf('initial guess was better: %g, (x,y)=(%g,%g), s=%g\n', ...
            mindd_guess - dist^2, x, y, s);
    tg = [xp(s_guess), yp(s_guess)];
    nor = [x - xs(s_guess), y-ys(s_guess)];
    figure(10);
    plot([x xs(s)], [y ys(s)], 'b.-');
    figure(11);
    plot(s, d2(s,x,y), 'b*');
    figure(12);
    s2 = s_guess
    plot(s2, g(s2,x,y), 'b*');
    plot([s2 s2+.27], g(s2,x,y) + [0  .27*gp(s2,x,y)], 'b-');
    grid on;
    figure(13);
    plot(s, gp(s,x,y), 'b*');
    keyboard
  end
end
end % function
