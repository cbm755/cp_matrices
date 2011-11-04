function [cpx,cpy,cpz, dist, bdy] = cpParamSurface(xx,yy,zz, paramf, paramf2nd, optfin, surfmesh, LB, UB, paramAdjust, How, DEBUG)
%CPPARAMSURFACE   CP representation of parameterised surfaces via optimization
%   Probably you don't want this directly, see e.g., cpMobiusStrip
%
%   [cpx,cpy,cpz dist, varargout] = cpParamSurface(p,q,r,
%                         x,xu,xv,xuu,xvv,xuv,
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
%    How: what technique to use, currently unsupported, might
%         eventually allow a selection from different implementations.
%
%
% TODO: this code is in a state of flux...  BEWARE

  if (nargin < 8)
    LB = [-inf -inf]
  end
  if (nargin < 9)
    UB = [inf inf]
  end
  if (nargin < 10)
    paramAdjust = @(my) my
  end
  if (nargin < 11)
    How = 0;
  end
  if (nargin < 12)
    DEBUG = 0;
  end

  if isempty(optfin)
    optf = @myfmincon_optf;
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

  [xp, yp, zp, up, vp] = surfmesh{:};

  for i = 1:nx
    xpt = [x1d(i); y1d(i); z1d(i)];
    p = xpt(1);
    q = xpt(2);
    r = xpt(3);

    %% initial guess
    A = tic;
    dds = (xp - xpt(1)).^2 + (yp - xpt(2)).^2 + (zp - xpt(3)).^2;
    [mindd_guess,I] = min(dds(:));
    s_initial_guess = [up(I); vp(I)];
    %[s_initial_guess, mindd_guess] = helper_initialguess(xpt, surfmesh);
    time_guess = toc(A);

    %
    out = [];

    if (How == 0)
      A = tic;
      [cp, dist1, bdy1] = helper_fmincon(xpt, optf, s_initial_guess, paramf, fmincon_opt, LB, UB, paramAdjust);
      out = [out toc(A)];
    end


    if (How == 1)
      % todo: in my tests this was slower, but should probaby fix it up anyway
      A = tic;
      [cp, dist1, bdy1] = helper_lsqnonlin(xpt, myvec, f, s_initial_guess, x, lsq_opt, LB, UB, paramAdjust);
      out = [out toc(A)];
    end


    if (How == 2)
      % todo: doesn't yet do boundaries, see comments
      A = tic;
      [cp, dist1, bdy1] = helper_newton(xpt, mindd_guess, s_initial_guess, paramf, paramf2nd);
      out = [out toc(A)];
    end

    %[time_guess out]

    cpx(i) = cp(1);
    cpy(i) = cp(2);
    cpz(i) = cp(3);
    dist(i) = dist1;
    bdy(i) = bdy1;
  end

  cpx = reshape(cpx, size(xx));
  cpy = reshape(cpy, size(xx));
  cpz = reshape(cpz, size(xx));
  dist = reshape(dist, size(xx));
  bdy = reshape(bdy, size(xx));

  %% helper functions placed inside: access local variables
  % careful: variables inside this aren't local either (!)
  function [d2, grad] = myfmincon_optf(uv, p)
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

end % end main function



%% helper functions



function [cp, dist, bdy] = helper_newton(xpt, mindd_guess, s_guess, paramf, paramf2nd)
% Newton's method
%
% Note: passing fcns to compute F and Jacobian was slow (or more likely
% creating those fcns for each xpt)
%
% TODO: need to deal properly with boundaries.  Either google a bit,
% or do a separate search of a parameterized boundary...
% Straightforward for mobius, just need to pass in some other
% functions, maybe in cell array that gets forwarded to the helper
% function

  DEBUG = 0;

  tol = 1e-14;
  maxn = 100;

  n = 1;
  outsideCounter = 0;
  s = s_guess;

  p = xpt(1);
  q = xpt(2);
  r = xpt(3);

  while (true)
    n = n + 1;
    %% Newton's method
    %J2 = Jfh(s);
    %norm(J)
    %max(max(J))
    %f2 = fh(s);
    %[x, xu, xv, xuu, xvv, xuv] = paramf(s(1), s(2));
    [x, xu, xv] = paramf(s(1), s(2));
    [xuu, xvv, xuv] = paramf2nd(s(1), s(2));

    f = [ 2*(x(1)-p) .* xu(1) + ...
          2*(x(2)-q) .* xu(2) + ...
          2*(x(3)-r) .* xu(3); ...
          2*(x(1)-p) .* xv(1) + ...
          2*(x(2)-q) .* xv(2) + ...
          2*(x(3)-r) .* xv(3) ];

    f1u = 2*(x(1)-p) .* xuu(1)  +  2*xu(1).^2 + ...
          2*(x(2)-q) .* xuu(2)  +  2*xu(2).^2 + ...
          2*(x(3)-r) .* xuu(3)  +  2*xu(3).^2;

    f2v = 2*(x(1)-p) .* xvv(1)  +  2*xv(1).^2 + ...
          2*(x(2)-q) .* xvv(2)  +  2*xv(2).^2 + ...
          2*(x(3)-r) .* xvv(3)  +  2*xv(3).^2;

    f1v = 2*(x(1)-p) .* xuv(1)  +  2*xu(1) .* xv(1) + ...
          2*(x(2)-q) .* xuv(2)  +  2*xu(2) .* xv(2) + ...
          2*(x(3)-r) .* xuv(3)  +  2*xu(3) .* xv(3);


    gamma = 10;
    if (s(2) > 1)
      f(2) = f(2) + gamma*(s(2)-1);
      f2v = f2v + gamma;
    elseif (s(2) < -1)
      f(2) = f(2) + gamma*(s(2)+1);
      f2v = f2v + gamma;
    end
    %penalty = (v<-1) .* gamma/2*((v+1).^2)  +  (v>1) .* gamma/2*((v-1).^2);
    % add that to the d2 function?  xx here is v

    J = [f1u f1v; f1v f2v];
    %[f-f2  J2-J]

    % J (snew-s) = -f

    snew  = s + J \ -f;

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
    % smooth at this point
    snew = paramAdjust(snew);

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
      if (abs(snew(2)) > 1)
        fprintf('converged: n=%d, s=(%f,%f)\n',n, snew(1),snew(2));
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

  %% TODO: not right
  if (abs(s(2)) > 1)
    bdy = true; % at least would need to search the bdy and I'm still
             % not convinced.
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
end % main function



function [F, J] = mylsqfun(s, p)
  [x, xu, xv, xuu, xvv, xuv] = mobius_parm(s(1), s(2), 1, 0.35);
  F = [x(1) - p(1); x(2) - p(2); x(3) - p(3)];
  if nargout > 1   % Two output arguments
    J = [xu(1)  xv(1); xu(2)  xv(2); xu(3)  xv(3)];
  end
end

function [cp, dist, bdy] = helper_lsqnonlin(xpt, d2, grad, initial_guess, paramf, opt, LB, UB, paramAdjust)

  t0 = initial_guess;

  %[res, fval, fvals, flg, output,lambda]...
  %    = lsqnonlin(d2, t0, [-inf; -1], [inf; 1], opt);
  [res, fval, fvals, flg, output,lambda]...
      = lsqnonlin(@(s) mylsqfun(s,xpt), t0, LB, UB, opt);

  if abs(abs(res(2)) - 1) < 1e-12
    %disp('** found boundary **');
    bdy = true;
  else
    %disp('** not boundary **');
    bdy = false;
  end

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
    warning('final search: possibly noncoverged CP search');
    keyboard
  end

  newres = paramAdjust(res);
  cp = paramf(newres(1),newres(2));
  dd = fval;
  % for lsqnonlin, need squared
  if (abs(dd - norm(cp - xpt,2)^2) > 1e-15)
    dd
    norm(cp - xpt,2)^2
    dd - norm(cp - xpt,2)^2
    warning('actual distance doesn''t match opt result')
    keyboard
  end

  dist = sqrt(dd);

  if (MakePlots)
    plot3([xpt(1) cp(1)], [xpt(2) cp(2)], [xpt(3) cp(3)], 'bo-')
    drawnow()
  end
  %[time_guess time_opt output.iterations]

end




function [cp, dist, bdy] = helper_fmincon(xpt, f, initial_guess, paramf, opt, LB, UB, paramAdjust)
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


  %[res(1) res(2) fval flg lambda.lower(2) lambda.upper(2)]

  % check lagrange mult
  if (lambda.lower(2) ~= 0)
    bdy11 = 1;
  else
    bdy11 = 0;
  end
  if (lambda.upper(2) ~= 0)
    bdy12 = 1;
  else
    bdy12 = 0;
  end
  bdy = bdy11 | bdy12;

  % double-check the lambda stuff: TODO: remove later
  if abs(abs(res(2)) - 1) < 1e-12
    if (~bdy)
      disp('** mismatch, lagrange mult didn''t find bdy? **');
      keyboard
    end
  else
    if (bdy)
      disp('** mismatch, lagrange mult found bdy? **');
      keyboard
    end
  end

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

  newres = paramAdjust(res);
  cp = paramf(newres(1),newres(2));
  dist = fval;
  if (abs(dist - norm(cp - xpt,2)^2) > 2e-15)
    dist
    norm(cp - xpt,2)^2
    dist - norm(cp - xpt,2)^2
    warning('actual distance doesn''t match opt result')
    keyboard
  end

  if (MakePlots)
    plot3([xpt(1) cp(1)], [xpt(2) cp(2)], [xpt(3) cp(3)], 'bo-')
    drawnow()
  end
  %[time_guess time_opt output.iterations]

end
