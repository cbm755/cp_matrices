function [cpx, cpy, dist, varargout] = cpParamCurveClosed_oldloop(x,y,xs,ys,xp,yp,xpp,ypp,endpt,DEBUG)
% matlab version of function in python code ParamCurve
% Find closest points of parameterised curves using Newton's method

% x,y : point to find closest point for
% xs,ys : parameterised curve
% xp,yp : first derivative of curve function
% xpp,ypp : second derivative of curve function
% endpt = [pt1 pt2]: endvalues of parameter (e.g., [0 2*pi])
%
% TODO: this code may produce ominous sounding messages that are
% actually harmless.  Converse may also be the case :(

% TODO: currently in 2D, generalise to curves in 3D

%%

if (nargin < 10)
  DEBUG = 0;
end

% distance squared from (x,y) to curve
d2 = @(t,x,y) (xs(t) - x).^2 + (ys(t) - y).^2;

% derivative of dist wrt t
g = @(t,x,y) 2*(xs(t) - x).*xp(t) + 2*(ys(t) - y).*yp(t);

% second derivative:
gp = @(t,x,y) 2*xp(t).*xp(t) + 2*(xs(t) - x).*xpp(t) + 2*yp(t).*yp(t) + 2*(ys(t) - y).*ypp(t);

ss = linspace(endpt(1), endpt(2), 500);
dd = d2(ss, x, y);
[mindd_guess, i] = min(dd);
s_guess = ss(i);

if (DEBUG >= 1)
  figure(10); clf; hold on;
  tt = linspace(s_guess-.2,s_guess+.2,5000);
  plot(xs(tt),ys(tt),'k-');
  %plot(x,y,'k*');
  plot([x xs(s_guess)],[y ys(s_guess)],'r-*');
  axis equal

  figure(11); clf; hold on;
  plot(tt, d2(tt,x,y), 'k-');
  plot(ss,dd,'k.');
  plot(s_guess, d2(s_guess,x,y), 'r*');

  figure(12); clf; hold on;
  plot(tt, g(tt,x,y), 'k-');
  plot(s_guess, 0, 'r*');

  figure(13); clf; hold on;
  plot(tt, gp(tt,x,y), 'k-');
  plot(s_guess, gp(s_guess,x,y), 'r*');

  %plot(tt, d2(tt), 'k-');
  %plot([x xs(s)], [y ys(s)],'b-o');
end

%%
% Newton's method

tol = 1e-14;
endpttol = 0.5;
maxn = 100;

n = 1;
outsideCounter = 0;
s = s_guess;

while (true)
  n = n + 1;
  %% Newton's method
  if (abs(gp(s,x,y)) > tol)
    snew = s - g(s,x,y) / gp(s,x,y);
  else
    % second deriv is zero
    snew = s;
  end

  %% Force parameter to be periodic
  % TODO: Newton's method if probably not robust if the curve is not
  % smooth at this point
  snew2 = argPeriodic(snew-endpt(1), endpt(2)-endpt(1)) + endpt(1);
  %[s snew s-snew snew2 snew2-snew]
  snew = snew2;

  if (DEBUG >= 10)
    figure(10);
    plot([xs(s)], [ys(s)],'bo');
    drawnow();
  end

  if (n > maxn)
       fail = 1;
       fprintf('too many iterations: (n,snew,x,y) = %d, %f, %f, %f \n', n, snew, x, y);
       g(s,x,y)
       gp(s,x,y)
       warning('max iterations in Newton solve: CP is likely wrong!');
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

cpx = xs(s); cpy = ys(s);
dist = sqrt((cpx - x).^2 + (cpy - y).^2);

if (fail)
  varargout = {1};
else
  varargout = {0};
end

if (DEBUG <= 0)
  assert(mindd_guess+tol >= dist^2, ...
         'initial guess was better: %g, (x,y)=(%g,%g), s=%g\n', ...
         mindd_guess - dist^2, x, y, s);
else
  if ~(mindd_guess+tol >= dist^2)
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


function u = argPeriodic(v, P)
%ARGPERIODIC  Return the principal value the argument.
%   argPeriodic(v, P) return a value in [0, P] where P defaults to
%   2*pi if not provided.

  if (nargin == 1)
    P = 2*pi;
  end

  A = v/P;
  B = A - floor(A);
  u = B*P;
end
