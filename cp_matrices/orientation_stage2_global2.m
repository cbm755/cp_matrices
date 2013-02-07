function [inside,outside,unknown] = orientation_stage2_global2(xx,yy,zz,cpx,cpy,cpz,dist,dx,E,inside0,outside0,verbose)
%ORIENTATION_STAGE2_GLOBAL  Helper function for ORIENTATION_FROM_CP
%   Classify a few closest points as inside/outside supposing that
%   most are already known.
%
%   Uses a global calculation based on interpolating distances and
%   neighbouring CPs.
%
%   WARNING: not well tested

  %% Parameters
  nnbrs = 6;   % how many neighbours near the CP
  DROP_TO_KEYBOARD_OK = 1;
  % important that this stays same as in normals_from_cp.m
  TOO_CLOSE_TOL = 100*eps;

  nnc0 = nnz(~inside0 & ~outside0);
  unknown = ~inside0 & ~outside0;
  I = find(unknown);

  sdist = -1*inside0.*abs(dist) + outside0.*abs(dist);

  E1 = E(:,unknown);
  E2 = E(:,~unknown);
  b = E2 * sdist(~unknown);

  % We want to solve a integer programming problem: we know the
  % distances but don't know the signs.

  % use the least squares linear problem solution as an initial guess:
  x0 = E1 \ -b;

  % We try to solve the optimization problem with a penalty.  This
  % requires the optimization toolbox.
  obj = @(x) sum((E1*x + b).^2);
  %ineqcons = @(x) [];
  %eqcons = @(x) x.^2 - abs(dist(unknown)).^2;
  %cons = @(x) cellfun(@(f) f(x), {ineqcons, eqcons});

  cons = @(x) mycon(x, abs(dist(unknown)));
  opt = optimset('tolfun', 1e-13, 'tolcon', 1e-13);
  [x,fval,flag] = fmincon(obj, x0, [],[],[],[],[],[], cons, opt);
  %[fval  flag  min(abs(x)-dist(unknown))  max(abs(x)-dist(unknown))]

  x3 = x;
  for i=1:nnc0
    x2 = x;
    x2(i) = -1*x(i);
    if obj(x2) < fval
      disp(sprintf('%i, global lower %g - %g = %g', i, fval, obj(x2), ...
                   fval-obj(x2)))
      x3(i) = -1*x(i);
    end
  end

  fval3 = obj(x3);
  if (fval < fval3)
    error('global error not found');
  else
    x = x3;  % let's hope this is closer to global min
  end
  %sdist2 = sign(X) .* abs(dist(I));
  %[E1*sdist2 + b  E*dist]
  %E1 * sdist2 + b

  inside = inside0;
  outside = outside0;
  inside(unknown) = sign(x) <= 0;
  outside(unknown) = sign(x) > 0;
  sdist2 = -1*inside.*abs(dist)  +  outside.*abs(dist);
  unknown = [];

  return

  % TODO: split out to another function

  n1 = E*sdist2;
  % let's check each one using the real distance
  for i=1:nnc0
    sd1 = sdist2;
    sd2 = sd1;
    sd2(I(i)) = -1*sd1(I(i));
    v1=E*sd1;
    v2=E*sd2;
    % TODO: maybe we should check if they lower the global norm
    % instead of locally?
    if norm(v1) > norm(v2)
      disp(sprintf('%d lower global', i));
      disp([i  norm(v1)  norm(v2)  abs(v1(I(i)))  abs(v2(I(i)))  abs(v1(I(i)))-abs(v2(I(i))) ...
            x(i)  abs(sd1(I(i)))])
    end
    if abs(v1(I(i))) > abs(v2(I(i)))
      disp(sprintf('%d lower local', i));
      % swap
      %x(i) = -1*x(i);
      %disp([i norm(v1) norm(v2) norm(v1)-norm(v2) abs(v1(I(i))) ...
      %      abs(v2(I(i)))  abs(v1(I(i)))-abs(v2(I(i)))])
      disp([i norm(v1) norm(v2)  abs(v1(I(i)))  abs(v2(I(i)))  abs(v1(I(i)))-abs(v2(I(i))) ...
            x(i)  abs(sd1(I(i)))])
      %figure(11);
      %alpha(0.95);
      %plot3(xx(I(i)), yy(I(i)), zz(I(i)), 'k*');
      %[xp,yp,zp]=sphere(20);
      %xp = 0.1*dx*xp + xx(I(i));
      %yp = 0.1*dx*yp + yy(I(i));
      %zp = 0.1*dx*zp + zz(I(i));
      %surf(xp,yp,zp,4*ones(size(xp)));
      %shading flat
      %keyboard
    end
  end
  inside = inside0;
  outside = outside0;
  inside(unknown) = sign(x) <= 0;
  outside(unknown) = sign(x) > 0;
  sdist2 = -1*inside.*abs(dist)  +  outside.*abs(dist);

  n2 = E*sdist2;
  norm(n1) - norm(n2)
  unknown = [];

  keyboard
end


function [c,ceq] = mycon(x,a)
  c = [];
  ceq = x.^2 - a.^2;
end