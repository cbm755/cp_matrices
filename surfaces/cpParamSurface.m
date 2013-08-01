function [cpx,cpy,cpz, dist, bdy, uu, vv] = cpParamSurface(xx,yy,zz, paramf, paramf2nd, paramfEdge, optfin, surfmesh, LB, UB, paramAdjust, How, DEBUG)
%CPPARAMSURFACE   CP representation of parameterised surfaces via optimization
%   Currently for open surfaces with a single edge (e.g., mobius strip)
%   Probably you don't want this directly, see e.g., cpMobiusStrip
%
%   [cpx,cpy,cpz dist, varargout] = cpParamSurface(p,q,r,
%                         paramf, paramf2nd, edgeparamf.
%                         optf, surfmesh, LB, UB, paramAdjust,
%                         How, Debug)
%    p,q,r: points to find closest points of.
%
%    paramf: a function of the form [x,xu,xv] = paramf(u,v), which
%            evaluates the surface x(u,v) and xu and xv are the
%            partial w.r.t. u and v.
%
%    paramf2nd: similar to above but provides [xuu,xvv,xuv], the
%               2nd partials.
%
%    paramfEdge: a function which parameterizes an edge of an open
%                surface.  TODO: what should we do for a closed
%                surface?  []?
%                Currently there can only be one of these (the surface
%                can have a single edge, e.g., hemisphere or mobius
%                strip).
%
%    optf: override the detault optimization function.  pass [] to
%          use the default.
%
%    surfmesh: a mesh of the surface used to get a good initial guess,
%              stored in a cell array as {xp,yp,zp,up,vp}.  Should
%              have a good coverage of the surface, may 500 points?
%
%    LB, UB: lower/upper bounds for parameters.  If you have a
%            periodic parameter, see the next option.  TODO
%
%    Optional inputs:
%    paramAdjust: a function that is called after the optimization (or
%                 perhaps during in some cases, e.g., after each
%                 step).  CAREFUL WITH THIS, DON'T DO ARBITRARY THINGS
%                 WITH IT.
%
%    How: what technique to use, various implementations, see source.
%
%
% TODO: this code is in a state of flux...  BEWARE

  if (nargin < 9)
    LB = [-inf -inf]
  end
  if (nargin < 10)
    UB = [inf inf]
  end
  if (nargin < 11)
    paramAdjust = [];
  end
  if (nargin < 12)
    How = 0;
  end
  if (nargin < 13)
    DEBUG = 0;
  end

  if isempty(optfin)
    if How == 0
      optf = @my_fmincon_optf;
    elseif How == 1
      optf = @my_lsq_optf;
    else
      % unneeded?
    end
  else
    optf = optfin;
  end

  % seems like less accurate in the x direction is normal
  opt = optimset('tolfun', 1e-14, 'tolx', 1e-12, 'tolcon', 1e-13, ...
                 'maxfunevals', 100000, 'maxiter', 100000);
  opt = optimset(opt, 'Display', 'off');
  opt = optimset(opt, 'largescale', 'off');
  opt = optimset(opt, 'algorithm', 'active-set');
  opt = optimset(opt, 'GradObj','on');
  fmincon_opt = opt;

  % seems like less accurate in the x direction is normal
  opt = optimset('tolfun', 1e-14, 'tolx', 1e-12, 'tolcon', 1e-13, ...
                 'maxfunevals', 100000, 'maxiter', 100000);
  opt = optimset(opt, 'Display', 'off');
  opt = optimset(opt, 'Jacobian', 'on');
  lsq_opt = opt;



  %% loop over the points and find the closest point for each
  x1d = xx(:); y1d = yy(:); z1d = zz(:);
  nx = length(x1d);     % number of points

  % allocate space
  cpx = zeros(nx,1);
  cpy = zeros(nx,1);
  cpz = zeros(nx,1);
  dist = zeros(nx,1);
  bdy = zeros(nx,1);
  uu = zeros(nx,1);
  vv = zeros(nx,1);

  [xp, yp, zp, up, vp] = surfmesh{:};

  fprintf('cpParamSurface: starting to process %d points\n', nx);
  for i = 1:nx
    xpt = [x1d(i); y1d(i); z1d(i)];
    p = xpt(1);
    q = xpt(2);
    r = xpt(3);

    %% initial guess
    time_guess = cputime();
    % vectorized:
    dds = (xp - xpt(1)).^2 + (yp - xpt(2)).^2 + (zp - xpt(3)).^2;
    [mindd_guess,I] = min(dds(:));
    s_initial_guess = [up(I); vp(I)];
    %[s_initial_guess, mindd_guess] = helper_initialguess(xpt, surfmesh);
    time_guess = cputime() - time_guess;

    out = [];

    if (How == 0)
      opt_time = cputime();
      [cp, dist1, bdy1, s] = helper_fmincon(xpt, optf, s_initial_guess, paramf, fmincon_opt, LB, UB, paramAdjust);
      u1 = s(1);  v1 = s(2);
      out = [out  cputime() - opt_time];
    end


    if (How == 1)
      % in my tests this was much slower
      opt_time = cputime();
      [cp, dist1, bdy1, s] = helper_lsq(xpt, optf, s_initial_guess, paramf, lsq_opt, LB, UB, paramAdjust);
      u1 = s(1);  v1 = s(2);
      out = [out  cputime() - opt_time];
    end


    if (How == 2)
      % fastest by an order of magnitude, maybe less reliable, lots
      % of parameters to tune (!)
      opt_time = cputime();
      [cp, dist1, bdy1, s] = helper_newton(xpt, mindd_guess, s_initial_guess, paramf, paramf2nd, LB, UB, paramAdjust);
      u1 = s(1);  v1 = s(2);
      if (bdy1 == 1)
        [cpx2,cpy2,cpz2,dist2,s2] = cpParam3DCurveClosed(xpt(1), xpt(2), xpt(3), paramfEdge, [0 4*pi]);
        cp = [cpx2; cpy2; cpz2];
        dist1 = dist2;
        u1 = s2;
        v1 = 1;   % TODO: hardcoded for mobius
      end
      out = [out  cputime() - opt_time];
    end

    %disp([time_guess out])

    cpx(i) = cp(1);
    cpy(i) = cp(2);
    cpz(i) = cp(3);
    dist(i) = dist1;
    bdy(i) = bdy1;
    uu(i) = u1;
    vv(i) = v1;
  end

  cpx = reshape(cpx, size(xx));
  cpy = reshape(cpy, size(xx));
  cpz = reshape(cpz, size(xx));
  dist = reshape(dist, size(xx));
  bdy = reshape(bdy, size(xx));
  uu = reshape(uu, size(xx));
  vv = reshape(vv, size(xx));

  %% helper functions placed inside: access local variables
  % careful: variables inside this aren't local either (!)
  function [d2, grad] = my_fmincon_optf(uv, p)
    % todo, nargout here and below, faster when no grad? (probably not)
    [tx, txu, txv] = paramf(uv(1), uv(2));

    d2 = sum( (tx - p).^2 );
    %d2 = (x(1) - p(1))^2 + ...
    %     (x(2) - p(2))^2 + ...
    %     (x(3) - p(3))^2;

    grad = [ 2*(tx - p)' * (txu);
             2*(tx - p)' * (txv) ];

    % TODO:
    %grad = [ 2*(x(1)-p(1)) * xu(1) + ...
    %         2*(x(2)-p(2)) * xu(2) + ...
    %         2*(x(3)-p(3)) * xu(3); ...
    %         2*(x(1)-p(1)) * xv(1) + ...
    %         2*(x(2)-p(2)) * xv(2) + ...
    %         2*(x(3)-p(3)) * xv(3) ];
  end

  function [F, J] = my_lsq_optf(uv, p)
    [tx, txu, txv] = paramf(uv(1), uv(2));
    F = tx - p;
    %if nargout > 1   % Two output arguments
    J = [txu(1)  txv(1); txu(2)  txv(2); txu(3)  txv(3)];
    %end
  end

end  % end main function



%% helper functions


function [cp, dist, bdy, s] = helper_newton(xpt, mindd_guess, s_guess, paramf, paramf2nd, LB, UB, paramAdjust)
% Newton's method
%
% Note: passing fcns to compute F and Jacobian was slow (or more likely
% creating those fcns for each xpt)
%
% This doens't deal with boundaries---just puts a penalty to try to
% stop it from converging too far outside.  bdy will be set to 1 if
% it finishes on or outside the boundary.  You could then do a
% search on the boundary curve.

  DEBUG = 1;

  tol = 1e-14;
  maxn = 100;

  n = 0;

  %s_guess = [3.51; -1]
  s = s_guess;

  while (true)
    n = n + 1;
    %J2 = Jfh(s);
    %f2 = fh(s);
    %[x, xu, xv, xuu, xvv, xuv] = paramf(s(1), s(2));
    [x, xu, xv] = paramf(s(1), s(2));
    [xuu, xvv, xuv] = paramf2nd(s(1), s(2));

    f = [ 2*(x-xpt)' * xu; ...
          2*(x-xpt)' * xv ];

    f1u = 2*(x-xpt)' * xuu + 2*xu' * xu;
    f2v = 2*(x-xpt)' * xvv + 2*xv' * xv;
    f1v = 2*(x-xpt)' * xuv + 2*xu' * xv;

    % add a penalty to distance for outside the boundaries
    % TODO: lots of choice here and somewhat hardcoded for mobius
    % strip (e.g., boundary is only in v).
    gamma = 100;
    pow = 3;
    % (x(u,v) - xpt)^2 + gamma*(v-1)^3
    if (s(2) > UB(2))
      f(2) = f(2) + gamma*(s(2)-UB(2))^(pow-1);
      f2v = f2v + (pow-1)*gamma*(s(2)-UB(2))^(pow-2);
    elseif (s(2) < LB(2))
      f(2) = f(2) - gamma*(s(2)-LB(2))^(pow-1);
      f2v = f2v - (pow-1)*gamma*(s(2)-LB(2))^(pow-2);
    end
    %penalty = (v<-1) .* gamma/2*((v+1).^2)  +  (v>1) .* gamma/2*((v-1).^2);
    % add that to the d2 function?  xx here is v

    J = [f1u f1v; f1v f2v];

    %% Newton's method
    % we want to minimize d^2 by solving (f := grad d^2) = 0.  At
    % each step we want solve solve
    %    J (snew-s) = -f
    % that is
    %    snew  = s + (J \ -f);
    % (provided that moves in a decent direction)

    % first compute the change
    change = (J \ -f);

    % now make sure its a decent direction of d^2
    innerprod = -f' * change;
    if (innerprod < 0)
      % negative (/w small tolerance), so take a small step in the
      % gradient direction (should really do a line search)
      if (DEBUG > 1)
        fprintf('iter n=%d, taking gradient decent direction, x=[%f,%f,%f]\n', n, xpt(1), xpt(2), xpt(3));
      end
      thechange = -f * 1/norm(f) * 0.05;  % TODO: do better!
    else
      % we have a decent direction, now limit how far we can go
      changetol = 0.1;  % TODO: another parameter!
      normch = norm(change);
      if (normch > changetol)
        alpha = changetol*1/normch;
        if (DEBUG>=2)
          fprintf('constraining stepsize: n=%d, s=[%g,%g], normchange=%g, change=[%g,%g], x=[%f,%f,%f]\n',n, s(1), s(2), normch, change(1), change(2), xpt(1), xpt(2), xpt(3));
        end
      else
        alpha = 1;
      end
      thechange = alpha*change;
    end
    if (DEBUG>=2)
      fprintf('iter: n=%d, s=[%g,%g], thechange=[%g,%g], x=[%f,%f,%f]\n',n, s(1), s(2), thechange(1), thechange(2), xpt(1), xpt(2), xpt(3));
    end
    snew = s + thechange;
    if (DEBUG>=2)
      [xtemp] = paramf(snew(1), snew(2));
      dd = sum( (xtemp-xpt).^2 );
      [xtemp] = paramf(s_guess(1), s_guess(2));
      dd0 = sum( (xtemp-xpt).^2 );
      dt = 0.01;
      [xtemp] = paramf(s_guess(1)+dt, s_guess(2));  dd1 = sum( (xtemp-xpt).^2 );
      [xtemp] = paramf(s_guess(1)-dt, s_guess(2));  dd2 = sum( (xtemp-xpt).^2 );
      [xtemp] = paramf(s_guess(1), s_guess(2)+dt);  dd3 = sum( (xtemp-xpt).^2 );
      [xtemp] = paramf(s_guess(1), s_guess(2)-dt);  dd4 = sum( (xtemp-xpt).^2 );

      fprintf('  , w/ new change, dd=%d, dd0=%g, ddnbrs=[%g,%g,%g,%g], grad dd=[%g,%g]\n', dd, dd0, dd1,dd2,dd3,dd4, (dd1-dd2)/(2*dt), (dd3-dd4)/(2*dt));
      end

    %if (snew(2) > 1.5)
    %  change = 0.5*change;
    %  snew = s + change;
    %end


    %keyboard

    %if (abs(gp(s,x,y)) > tol)
    %TODO
    %gp(s,x,y) * (snew - s) =  - g(s,x,y);
    %else
    % second deriv is zero
    %snew = s;
    %end

    %% Do any adjustments to the parameter values
    % e.g, Force parameter to be periodic
    % TODO: Newton's method if probably not robust if the curve is not
    % smooth at this point, in general be careful with this.
    if (~isempty(paramAdjust))
      snew = paramAdjust(snew);
    end

  if (DEBUG >= 10)
    cp = paramf(snew(1), snew(2));
    fprintf('iter: n=%d, snew=(%f,%f), s0=(%f,%f), x=(%f,%f,%f)\n',n, snew(1), snew(2), s_guess(1), s_guess(2), xpt(1), xpt(2), xpt(3))
    %figure(1);
    set(0, 'CurrentFigure', 1);
    plot3(cp(1),cp(2),cp(3),'bo');
    axis equal
    drawnow();
    pause
  end

    %converged: n=52, s=(3.133150,1.011001), x=(-0.343126,0.148547,0.514520)
    %too many iterations: (n,snew,x) = 201, (8.958740,0.999768), (-0.333417,-0.101126,-0.078792)

    if (n > maxn)
      fail = 1;
      sdiffnorm = norm(s-s_guess);
      fprintf('too many iters: n=%d, sdiff=%f, s=(%f,%f), s0=(%f,%f), x=(%f,%f,%f)\n',n, sdiffnorm, snew(1), snew(2), s_guess(1), s_guess(2), xpt(1), xpt(2), xpt(3));
      %fprintf('too many iterations: (n,snew,x) = %d, (%f,%f), (%f,%f,%f)\n', n, snew(1),snew(2), xpt(1), xpt(2), xpt(3));
      %fprintf('  (s0) = (%f,%f)\n', s_guess(1), s_guess(2));
      %g(s,x,y)
      %gp(s,x,y)
      warning('max iterations in Newton solve: CP is likely wrong!');
      keyboard
      break;
    elseif (abs(s-snew) < tol)
      fail = 0;
      if ( ((abs(snew(2)) > 1) && (abs(snew(2)) > 1.2))  ||  ...
           (n > 20) )
        sdiffnorm = norm(s-s_guess);
        fprintf('converged: n=%d, sdiff=%f, s=(%f,%f), s0=(%f,%f), x=(%f,%f,%f)\n',n, sdiffnorm, snew(1), snew(2), s_guess(1), s_guess(2), xpt(1), xpt(2), xpt(3));
        %keyboard
      end
      %fprintf('converged, n=%d\n', n);
      break;
    end
    % update
    s = snew;
  end

  s = snew;
  %cp = [xparam{1}(s(1),s(2));  xparam{2}(s(1),s(2)); xparam{3}(s(1),s(2))];
  cp = paramf(s(1), s(2));
  %dist = sqrt((cpx - p).^2 + (cpy - q).^2 + (cpz - r).^2);
  dist = norm(xpt - cp, 2);


  % somewhat hardcoded for mobius (no s1 here)
  if ( (s(2) < LB(2)) || (s(2) > UB(2)) )
    bdy = true;
    % need to search the boundary (we haven't found it properly yet)
  else
    bdy = false;
  end

  if (fail)
    varargout = {1};
  else
    varargout = {0};
  end

  if (DEBUG <= 0)
    %assert(mindd_guess+tol >= dist^2, ...
    %       'initial guess was better: %g, (x,y)=(%g,%g,%g), s=%g\n', ...
    %       mindd_guess - dist^2, xpt(1),xpt(2),xpt(3), s);
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
end  % newton function





function [cp, dist, bdy, res] = helper_lsq(xpt, f, initial_guess, paramf, opt, LB, UB, paramAdjust)

  t0 = initial_guess;

  %[res, fval, fvals, flg, output,lambda]...
  %    = lsqnonlin(d2, t0, [-inf; -1], [inf; 1], opt);
  [res, fval, fvals, flg, output,lambda]...
      = lsqnonlin(@(s) f(s,xpt), t0, LB, UB, opt);

  bdy = any([lambda.lower; lambda.upper] ~= 0);
  % TODO: lambda more reliable?
  %if abs(abs(res(2)) - 1) < 1e-12
  %  bdy = true;
  %else
  %  bdy = false;
  %end
  %[bdy lambda.lower' lambda.upper']

  if (flg < 1)
    xpt, res, t0, fval, fvals, flg, output
    lambda
    warning('lsq search: possibly noncoverged CP search');
    keyboard
  end

  if (~isempty(paramAdjust))
    newres = paramAdjust(res);
  else
    newres = res;
  end

  cp = paramf(newres(1),newres(2));
  dd = fval;
  % for lsqnonlin, need squared dist
  if (abs(dd - sum((cp - xpt).^2)) > 1e-15)
    dd
    sum((cp - xpt).^2)
    dd - sum((cp - xpt).^2)
    warning('actual distance doesn''t match opt result')
    keyboard
  end

  res = newres;
  dist = sqrt(dd);
end  % helper_lsq function




function [cp, dist, bdy, res] = helper_fmincon(xpt, f, initial_guess, paramf, opt, LB, UB, paramAdjust)
  MakePlots = 0;

  if (MakePlots)
    fn = 99;
    %figure(fn); clf;
    set(0, 'CurrentFigure', fn); clf;
    [x,y,z] = paramMobiusStrip(32,Rad,Thick);
    surf(x,y,z);
    xlabel('x'); ylabel('y'); zlabel('z');
    hold on;
    axis equal
  end

  if (MakePlots)
    cp_guess = [xp(I) yp(I) zp(I)];
    plot3([xpt(1) cp_guess(1)], [xpt(2) cp_guess(2)], ...
          [xpt(3) cp_guess(3)], 'ro--')
  end

  t0 = initial_guess;

  %[res, fval, flg, output,lambda,jacobian]...
  %    = fmincon(f, t0, ...
  %              [],[], [],[], [-inf; -1], [inf; 1], [],  opt);
  [res, fval, flg, output,lambda,jacobian]...
      = fmincon(@(s) f(s,xpt), t0, ...
                [],[], [],[], LB, UB, [],  opt);

  bdy = any([lambda.lower; lambda.upper] ~= 0);

  %[res(1) res(2) fval flg lambda.lower(2) lambda.upper(2)]
  % check lagrange mult
  %if (lambda.lower(2) ~= 0)
  %  bdy11 = 1;
  %else
  %  bdy11 = 0;
  %end
  %if (lambda.upper(2) ~= 0)
  %  bdy12 = 1;
  %else
  %  bdy12 = 0;
  %end
  %bdy = bdy11 | bdy12;

  % double-check the lambda stuff: TODO: remove later
  %if abs(abs(res(2)) - 1) < 1e-12
  %  if (~bdy)
  %    disp('** mismatch, lagrange mult didn''t find bdy? **');
  %    keyboard
  %  end
  %else
  %  if (bdy)
  %    disp('** mismatch, lagrange mult found bdy? **');
  %    keyboard
  %  end
  %end

  if (flg < 1)
    xpt
    res
    t0
    fval
    fvals
    flg
    output
    lambda
    %grad
    warning('possibly noncoverged CP search');
    keyboard
  end

  if (~isempty(paramAdjust))
    newres = paramAdjust(res);
  else
    newres = res;
  end
  cp = paramf(newres(1),newres(2));
  dd = fval;
  if (abs(dd - sum((cp - xpt).^2)) > 2e-15)
    dd
    sum((cp - xpt).^2)
    dd - sum((cp - xpt).^2)
    warning('actual distance doesn''t match opt result')
    keyboard
  end

  dist = sqrt(dd);
  res = newres;

  if (MakePlots)
    plot3([xpt(1) cp(1)], [xpt(2) cp(2)], [xpt(3) cp(3)], 'bo-')
    drawnow()
  end
end  % helper_fmincon function




