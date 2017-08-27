function [cpx,cpy,dist,bdy] = cpCutHole2d(x,y,cpf,hole_cen,arclen_target,tol,dx,bw,given_cp_data)
%cpCutHole2d  Cut holes in surfaces described by cpfun
%
% x,y coordinate matrices of any same dimension
%
% hole_cen = [hx hy] gives the center of the hole.  If this is not on
% the surface given by cpf, it will be projected onto the surface (by
% calling cpf).
%
% Multiple holes are supported: pass each center as a row of "hole_cen"
% and pass a vector to "arclen_target".
%
% See cpCutHole for caveats and additional documentation.

  if (nargin >= 9)
    cpx = given_cp_data{1};
    cpy = given_cp_data{2};
  else
    [cpx, cpy] = cpf(x, y);
    %[cpx, cpy, cpz, tilde, tilde] = cpf(x, y, z);
  end
  bdy = zeros(size(cpx));
  [num_holes, temp] = size(hole_cen);
  assert (temp == 2)

  %code below assumes cpx, cpy are long vectors
  savesz = size(cpx);
  cpx = cpx(:);
  cpy = cpy(:);
  bdy = bdy(:);

  for k=1:num_holes
    xh = hole_cen(k, 1)
    yh = hole_cen(k, 2)
    % project onto surface
    [xh yh] = cpf(xh, yh)

    % look for a sphere radius which gives approx the target arclength
    R_bar = arclen_target(k)/2
    RR = linspace(0.7*R_bar,1.02*R_bar, 100);
    for i = 1:length(RR)
      I = ((cpx-xh).^2 + (cpy-yh).^2 < RR(i)^2);
      assert(nnz(I) > 0, 'empty circle/surface intersection!')
      SASA(i) = sum(I(:))*(dx)/(2*bw);
    end
    %SASA./arclen_target(k)
    [tilde,i] = min(abs(SASA./arclen_target(k)-1))
    R_H = RR(i)
    assert (i > 1 && i < length(RR), 'should be in the middle')
    assert (tilde < 0.1)

    % index of all points whose cp are within the sphere
    I = ((cpx-xh).^2 + (cpy-yh).^2 < R_H^2);

    %% TODO: temp plotting
    figure(10);
    porcupine_plot2d(x, y, cpx, cpy);
    plot(xh, yh, 'r*', 'linewidth', 5);
    plot(cpx(I), cpy(I), 'r.');
    plot(x(I), y(I), 'bo');
    t = exp(1i*linspace(0,2*pi,64));
    plot(xh + R_bar*real(t), yh + R_bar*imag(t), 'c:');
    plot(xh + R_H*real(t), yh + R_H*imag(t), 'r-');


    %% fixed-point iteration
    xl = cpx(I);
    yl = cpy(I);

    sdist = 1;

    tol2 = 0.01;

    while (abs(sdist) > tol | abs(dist1) > tol)
      % project onto sphere "hole"
      [xl2, yl2, sdist] = cpCircle(xl, yl, R_H, [xh yh]);

      % project onto surface
      [xll, yll, dist1] = cpf(xl2, yl2);

      % avoid projecting back and forth infinite times
      dist2=sqrt((xl-xll).^2 + (yl-yll).^2);
      dist3=sqrt((xl-xl2).^2 + (yl-yl2).^2);
      if  dist2  < tol2 & dist3 > R_H/2
        xll = xll+10^(-5)*rand; % add perturbation
        yll = yll+10^(-5)*rand;
      end

      xl=xll;
      yl=yll;
    end

    cpx(I) = xl;
    cpy(I) = yl;

    bdy = bdy | I;
end


  cpx = reshape(cpx, savesz);
  cpy = reshape(cpy, savesz);
  bdy = reshape(bdy, savesz);
  dist = sqrt((cpx-x).^2 + (cpy-y).^2);
end
