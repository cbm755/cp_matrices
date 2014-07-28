function [cpx, cpy, dist, bdy, varargout] = cpParamCurveOpen(x,y,xs,ys,xp,yp,xpp,ypp,endpt,DEBUG)
% matlab version of function in python code ParamCurve
% Find closest points of parameterised curves using Newton's method

% x,y : point to find closest point for
% xs,ys : parameterised curve
% xp,yp : first derivative of curve function
% xpp,ypp : second derivative of curve function
% endpt = [pt1 pt2]: endpts for curve with boundaries / endvalues of parameter
%

% TODO: currently in 2D, generalise to curve in 3D

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


warning('This code doesn''t work (yet), see _oldloop version and cpSpline2D.m');


%% Initial guess
% use M equispaced samples in the parameter space
M = 100
tic
s_guess = zeros(size(x));
mindd_guess = zeros(size(x));
ss = linspace(endpt(1), endpt(2), M);

% split the point list into shorter chunks, for each of these we
% will process M points above, the matrix will be quite large.
totalN = length(x);
N = 500
numchunks = ceil(totalN/N)

SS = repmat(ss, N, 1);

for i=1:numchunks
  if i < numchunks
    I = (1+(i-1)*N):i*N;
  else
    I = (1+(i-1)*N):totalN;
    NN = length(I);
    SS = repmat(ss, NN, 1);
  end
  X = repmat(x(I), 1, M);
  Y = repmat(y(I), 1, M);
  dd = d2(SS, X, Y);
  [temp, ii] = min(dd, [], 2);
  s_guess(I) = ss(ii)';
  mindd_guess(I) = temp;
end
toc


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


tic
totalN = length(x);
N = 1000;
numchunks = ceil(totalN/N)

cpx = zeros(size(x));
cpy = zeros(size(x));
dist = zeros(size(x));
bdy = zeros(size(x));
gfail = zeros(size(x));

for i=1:numchunks
  if i < numchunks
    I = (1+(i-1)*N):i*N;
  else
    I = (1+(i-1)*N):totalN;
  end
  [cpx(I) cpy(I) dist(I) bdy(I) fail] = newton(g, gp, x(I), y(I), xs, ys, s_guess(I), mindd_guess(I), endpt, d2);
  if (fail)
    warning('one chunk failed');
    warning('TODO: could process in a loop here, or otherwise be more robust');
    gfail(I) = 1;
  else
    gfail(I) = 0;
  end
end
toc


end % function


function [cpx cpy dist fail] = newton(g, gp, x, y, xs, ys, g_guess, mindd_guess, endpt, d2)
%%
% Newton's method

DEBUG = 0;

tol = 1e-14;
endpttol = 0.5;
maxn = 100;

n = 1;
outsideCounter = 0;
s = s_guess;

Outsiders = zeros(size(x));

while (true)
  n = n + 1;

  % if the second-deriv of some components is close to zero, don't
  % update s for those ones.  Otherwise do a newton step.
  snew = s;
  I = abs(gp(s,x,y)) > tol;
  snew(I) = s(I) - g(s(I),x(I),y(I)) ./ gp(s(I),x(I),y(I));
  % don't update ones that have diverged
  I = (Outsiders > 5)
  snew(I) = s(I);

  OutSide = (snew - endpt(1) < -endpttol) | (snew - endpt(2) > endpttol);

  Outsiders = Outsiders + Outside;

  if ( (snew - endpt(1) < -endpttol) || (snew - endpt(2) > endpttol) )
    outsideCounter = outsideCounter + 1;
  end

   if (outsideCounter >= 5)
       fail = 1;
       if (DEBUG >= 1)
         fprintf(['*** probably diverging outside: (s_guess,n,s,snew)= ' ...
                  '%f, %d, %f, %f \n'], s_guess,n,s,snew);
       end
       break;
   elseif (n > maxn)
       fail = 1;
       if ( abs(s_guess - s) <= 0.1*(endpt(2)-endpt(1)) )
         fprintf('too many iterations: (n,s_guess,s,x,y) = %d,%g,%g,%g,%g\n', n, s_guess,snew, x, y);
         g(s,x,y)
         gp(s,x,y)
         fprintf(['(this can also indicate that the CP is an endpt, ' ...
                  'particularly if s is very different from ' ...
                  's_guess)\n']);
       end
       break;
   elseif (abs(s-snew) < tol)
       if ( (endpt(1) <= snew) && (snew <= endpt(2)) )
            fail = 0;
            %fprintf('converged: (cpx,cpy,s,n)= %f, %f, %f, %d \n',xs(snew),ys(snew),snew,n);
       else
           fail = 1;
           if (DEBUG >= 1)
             fprintf('converged outside interval (cpx,cpy,s,n)= %f, %f, %f, %d \n',xs(snew),ys(snew),snew,n);
           end
       end
       break;
   end
   %disp(n)
   % update
   s = snew;
end


if (fail)
   % not converged: use endpt value
   dpt1 = d2(endpt(1), x, y);
   dpt2 = d2(endpt(2), x, y);
   [dpt,i] = min([dpt1 dpt2]);
   dist = sqrt(dpt);
   cpx = xs(endpt(i)); cpy = ys(endpt(i));
   bdy = i;
else  % converged
   cpx = xs(s); cpy = ys(s);
   dist = sqrt((cpx - x).^2 + (cpy - y).^2);
   bdy = 0;
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
end
