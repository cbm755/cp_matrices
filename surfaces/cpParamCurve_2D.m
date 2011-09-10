function [cpx, cpy, dist, varargout] = cpParamCurve_2D(x,y,xs,ys,xp,yp,xpp,ypp,endpt1,endpt2,isclosed,DEBUG)
% matlab version of function in python code ParamCurve
% Find closest points of parameterised curves using Newton's method

% x,y : point to find closest point for
% xs,ys : parameterised curve
% xp,yp : first derivative of curve function
% xpp,ypp : second derivative of curve function
% endpt1,endpt2: endpts for curve with boundaries / endvalues of parameter
% isclosed: true for closed curves

% TODO: think about endpoints / curves without bndries

% TODO: currently in 2D, generalise to 3D

%%

if (nargin < 12)
  DEBUG = 0;
end

% distance squared from (x,y) to curve
d2 = @(t,x,y) (xs(t) - x).^2 + (ys(t) - y).^2;

% derivative of dist wrt t
g = @(t,x,y) 2*(xs(t) - x).*xp(t) + 2*(ys(t) - y).*yp(t);

% second derivative:
gp = @(t,x,y) 2*xp(t).*xp(t) + 2*(xs(t) - x).*xpp(t) + 2*yp(t).*yp(t) + 2*(ys(t) - y).*ypp(t);

ss = linspace(endpt1, endpt2, 500);
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
   
   % Newton's method
   if (abs(gp(s,x,y)) > tol)
       snew = s - g(s,x,y) / gp(s,x,y);
   else
     % second deriv is zero
     snew = s;
   end

   % if closed surface, force it to be periodic in the parameter
   if (isclosed)
   %if (sp.isPeriodic == 1)
   %disp('blah');
    snew2 = argPeriodic(snew-endpt1, endpt2-endpt1) + endpt1;
    %[s snew s-snew snew2 snew2-snew]
    snew = snew2;
    %else
    %tp = t;
  end

  if (DEBUG >= 10)
    figure(10);
    plot([xs(s)], [ys(s)],'bo');
    drawnow();
  end

  if (~isclosed && (snew - endpt1 < -endpttol) || (snew - endpt2 > endpttol))
    outsideCounter = outsideCounter + 1
  end

   if (outsideCounter >= 5)
       fail = 1;
       fprintf('*** probably diverging: (s_guess,n,s,snew)= %f, %d, %f, %f \n', s_guess,n,s,snew);
       break;
   elseif (n > maxn)
       fail = 1;
       fprintf('too many iterations: (n,snew,x,y) = %d, %f, %f, %f \n', n, snew, x, y);
       g(s,x,y)
       gp(s,x,y)
       break;
   elseif (abs(s-snew) < tol)
       if (isclosed || ((endpt1 <= snew) && (snew <= endpt2)))
            fail = 0;
            %fprintf('converged: (cpx,cpy,s,n)= %f, %f, %f, %d \n',xs(snew),ys(snew),snew,n);
       else
           fail = 1;
           fprintf('converged outside interval (cpx,cpy,s,n)= %f, %f, %f, %d \n',xs(snew),ys(snew),snew,n);
       end;
       break;
   end
   %disp(n)
   % update
   s = snew;
end


if (fail)
   % not converged: use endpt value
   dpt1 = d2(endpt1, x, y);
   dpt2 = d2(endpt2, x, y);
   [dpt,i] = min([dpt1 dpt2]);
   dist = sqrt(dpt);
   endpt = [endpt1, endpt2];
   cpx = xs(endpt(i)); cpy = ys(endpt(i));

   varargout = {1};

else   % converged

    cpx = xs(s); cpy = ys(s);
    dist = sqrt((cpx - x).^2 + (cpy - y).^2);

    varargout = {0};

end

if (DEBUG <= 0)
  assert(mindd_guess+tol >= dist^2, ...
         'initial guess was better: %g, (x,y)=(%g,%g), s=%g\n', ...
         mindd_guess - dist^2, x, y, s ;
else
  if ~(mindd_guess+tol >= dist^2)
    fprintf('initial guess was better: %g, (x,y)=(%g,%g), s=%g\n', ...
            mindd_guess - dist^2, x, y, s );
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
end

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
