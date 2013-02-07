function [inside,outside,unknown] = orientation_stage2_global1(xx,yy,zz,cpx,cpy,cpz,dist,dx,E,inside0,outside0,verbose)
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
  % distances but don't know the signs.  Instead we solve a least
  % squares problem for both distances and signs.  Then we'll extract
  % the sign.

  % TODO: there are some zero rows of E1 that should be excluded here
  X = E1 \ -b;

  inside = inside0;
  outside = outside0;
  inside(unknown) = sign(X) <= 0;   % call zero inside too
  outside(unknown) = sign(X) > 0;
  sdist2 = -1*inside.*abs(dist)  +  outside.*abs(dist);
  unknown = [];

  % TODO: we could stop here, or do some more tests
  return

  % TODO: we should split these other tests out to another function

  % We can check each one using the real distance.  In my limited
  % tests, swapping based on local gave the same results as the
  % other more local code.
  for i=1:nnc0
    sd1 = sdist2;
    sd2 = sd1;
    sd2(I(i)) = -1*sd1(I(i));
    v1=E*sd1;
    v2=E*sd2;
    % global decrease
    if norm(v1) > norm(v2)
      disp(sprintf('%d lower global', i));
      disp([i  norm(v1)  norm(v2)  abs(v1(I(i)))  abs(v2(I(i)))  abs(v1(I(i)))-abs(v2(I(i))) ...
            X(i)  abs(sd1(I(i)))])

    end
    % local decrease
    if abs(v1(I(i))) > abs(v2(I(i)))
      disp(sprintf('%d lower local', i));
      % TODO: here we swap to minimize the local error
      X(i) = -1*X(i);
      %disp([i norm(v1) norm(v2) norm(v1)-norm(v2) abs(v1(I(i))) ...
      %      abs(v2(I(i)))  abs(v1(I(i)))-abs(v2(I(i)))])
      disp([i norm(v1) norm(v2)  abs(v1(I(i)))  abs(v2(I(i)))  abs(v1(I(i)))-abs(v2(I(i))) ...
            X(i)  abs(sd1(I(i)))])
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
  inside(unknown) = sign(X) <= 0;
  outside(unknown) = sign(X) > 0;
  sdist2 = -1*inside.*abs(dist)  +  outside.*abs(dist);
  unknown = [];

