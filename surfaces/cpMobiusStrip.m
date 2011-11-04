function [cpx,cpy,cpz, dist, bdy, uu, vv] = cpMobiusStrip(x, y, z, Rad, Thick, cen)
%CPMOBIUSSTRIP  Closest point representation of a Mobius strip
%   [cpx,cpy,cpz, dist, bdy] = cpMobiusStrip(x, y, z)
%
%   [...] = cpMobiusStrip(x, y, z, Radius, Width, Center)

  % defaults
  if (nargin < 4)
    Rad = 1;
  end
  if (nargin < 5)
    Thick = 0.35;
  end
  if (nargin < 6)
    cen = [0 0 0];
  end

  %[xx, xxu, xxv, xxuu, xxvv, xxuv] = mobius_param_fcns(Rad, Thick);

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  paramf = @(u,v) mobius_parm(u, v, Rad, Thick);
  paramf2nd = @(u,v) mobius_parm_2ndpartials(u, v, Rad, Thick);
  optf = @mobius_distsqr_optf;
  surfmesh = {[] [] [] [] []};
  [surfmesh{:}] = paramMobiusStrip(100);
  paramAdjust = [];
  paramfEdge = @(u) mobius_edge_parm(u, 1, Rad, Thick);

  How = 2;
  [cpx,cpy,cpz,dist,bdy,uu,vv] = cpParamSurface(x,y,z, paramf, paramf2nd, paramfEdge, [], surfmesh, [-inf -1], [inf 1], paramAdjust, How);

  % TODO: should reset uu to [0,2*pi) here (and swap sign of v)



  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);

  function [dd, grad] = mobius_distsqr_optf(t, xpt)
    % return squared distance and its gradient, designed for "fmincon".
    u = t(1);
    v = t(2);
    %u = argPeriodic(u);

    [x,xu,xv] = mobius_parm(u, v, Rad, Thick);
    dd = sum( (x - xpt).^2 );
    %dd = norm(x - xpt, 2)^2;  % slow b/c sqrt!

    grad(1) = 2*(x - xpt)' * (xu);
    grad(2) = 2*(x - xpt)' * (xv);
  end

  function [F, J] = mobius_lsqfun(t, p)
    % return point difference vector (and Jacobian), designed for
    % "lsqnonlin".
    u = t(1);
    v = t(2);
    [x,xu,xv] = mobius_parm(u, v, Rad, Thick);
    F = x - p;
    %if (nargout > 1)   % Two output arguments
    J = [xu(1)  xv(1); xu(2)  xv(2); xu(3)  xv(3)];
    %end
end

end % main cpMobiusStrip



%% other helper functions

function uvout = mobius_param_adjust(uvin)
%Adjust parameters for Mobius Strip

  uvout = zeros(size(uvin));
  uvout(1) = argPeriodic(uvin(1), 4*pi);
  uvout(2) = uvin(2);
  % TODO: also possible to make 2*pi periodic and swap sign of v
  %cp = mobius_parm(newres(1),newres(2), Rad, Thick);
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




function [xx, xxu, xxv] = mobius_parm(u,v, R, T)
  %x = R*(1 + T*v .* cos(.5*u)) .* cos(u) - R*T/2;
  %y = R*(1 + T*v .* cos(.5*u)) .* sin(u);
  %z = R*T*v .* sin(0.5*u);
  t1 = R*(1 + T*v .* cos(.5*u));
  x = t1 .* cos(u) - R*T/2;
  y = t1 .* sin(u);
  z = R*T*v .* sin(0.5*u);
  xx = [x;y;z];
  xu = t1 .* -sin(u)   +   R*(T*v .* -sin(.5*u)*0.5) .* cos(u);
  yu = t1 .* cos(u)    +   R*(T*v .* -sin(.5*u)*0.5) .* sin(u);
  xv = R*T * cos(.5*u) .* cos(u);
  yv = R*T * cos(.5*u) .* sin(u);
  zu = R*T*v .* cos(0.5*u) * 0.5;
  zv = R*T * sin(0.5*u);
  xxu = [xu; yu; zu];
  xxv = [xv; yv; zv];
end

function [xxuu, xxvv, xxuv] = mobius_parm_2ndpartials(u,v, R, T)
  t1 = R*(1 + T*v .* cos(.5*u));
  %x = t1 .* cos(u) - R*T/2;
  %y = t1 .* sin(u);
  %z = R*T*v .* sin(0.5*u);
  %xx = [x;y;z];
  xuu = R*T*v.*sin(1/2*u).*sin(u) - 1/4*R*T*v.*cos(1/2*u).*cos(u) - R*cos(u).*(T*v.*cos(1/2*u) + 1);
  yuu = -R*sin(u)*(T*v.*cos(1/2*u) + 1) - 1/4*R*T*v.*cos(1/2*u).*sin(u) - R*T*v.*sin(1/2*u).*cos(u);
  zuu = -1/4*R*T*v.*sin(1/2*u);
  xvv = 0;
  yvv = 0;
  zvv = 0;
  xuv = -R*T*cos(1/2*u).*sin(u) - 1/2*R*T*sin(1/2*u).*cos(u);
  yuv = R*T*cos(1/2*u).*cos(u) - 1/2*R*T*sin(1/2*u).*sin(u);
  zuv = 1/2*R*T*cos(1/2*u);
  xxuu = [xuu; yuu; zuu];
  xxvv = [xvv; yvv; zvv];
  xxuv = [xuv; yuv; zuv];
end



function [xx, xxu, xxuu] = mobius_edge_parm(u, v, R, T)
  %v = 1;
  %x = R*(1 + T*v .* cos(.5*u)) .* cos(u) - R*T/2;
  %y = R*(1 + T*v .* cos(.5*u)) .* sin(u);
  %z = R*T*v .* sin(0.5*u);

  t1 = R*(1 + T*v .* cos(.5*u));
  x = t1 .* cos(u) - R*T/2;
  y = t1 .* sin(u);
  z = R*T*v .* sin(0.5*u);
  xx = [x;y;z];

  xu = t1 .* -sin(u)   +   R*(T*v .* -sin(.5*u)*0.5) .* cos(u);
  yu = t1 .* cos(u)    +   R*(T*v .* -sin(.5*u)*0.5) .* sin(u);
  zu = R*T*v .* cos(0.5*u) * 0.5;
  xxu = [xu; yu; zu];

  xuu = R*T*v.*sin(1/2*u).*sin(u) - 1/4*R*T*v.*cos(1/2*u).*cos(u) - R*cos(u).*(T*v.*cos(1/2*u) + 1);
  yuu = -R*sin(u).*(T*v.*cos(1/2*u) + 1) - 1/4*R*T*v.*cos(1/2*u).*sin(u) - R*T*v.*sin(1/2*u).*cos(u);
  zuu = -1/4*R*T*v*sin(1/2*u);
  xxuu = [xuu; yuu; zuu];
end
