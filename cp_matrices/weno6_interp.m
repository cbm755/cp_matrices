function w = weno6_interp(cp, f, x)
%WENO6_INTERP  nonlinear WENO interpolation in 2D/3D
%   TODO, following weno4_interp(), not implemented completely/

  warning('***not fully implemented, not tested at all***');

  [n1,dim] = size(x);
  if dim == 2
    Cache = weno6_interp2d_makecache(cp, f, x);
    w = weno6_interp2d(Cache, f);
  elseif dim == 3
    Cache = weno6_interp3d_makecache(cp, f, x);
    w = weno6_interp3d(Cache, f);
  else
    error('dim not implemented');
  end
end




function w = weno6_interp2d(Cache, f)
%WENO6_INTERP2D  nonlinear WENO interpolation 2D

  C = Cache.C;
  xi = Cache.xi;
  yi = Cache.yi;
  x = Cache.x;
  y = Cache.y;
  dx = Cache.dx;

  tic
  u = {};
  for j=1:6
    u{j} = helper1d( ...
        {C{1,j}*f, C{2,j}*f, C{3,j}*f, C{4,j}*f, C{5,j}*f, C{6,j}*f}, ...
        xi, dx, x);
  end
  w = helper1d({u{1}, u{2}, u{3}, u{4}, u{5}, u{6}}, yi, dx, y);
  toc
end




function w = weno6_interp3d(Cache, f)
%WENO6_INTERP3D  nonlinear WENO interpolation 3D

  C = Cache.C;
  xi = Cache.xi;
  yi = Cache.yi;
  zi = Cache.zi;
  x = Cache.x;
  y = Cache.y;
  z = Cache.z;
  dx = Cache.dx;

  tic
  u = {};
  v = {};
  for k=1:6
    for j=1:6
      u{j} = helper1d( ...
          {C{1,j,k}*f, C{2,j,k}*f, C{3,j,k}*f, C{4,j,k}*f, C{5,j,k}*f, C{6,j,k}*f}, ...
          xi, dx, x);
    end
    v{k} = helper1d({u{1}, u{2}, u{3}, u{4}, u{5}, u{6}}, yi, dx, y);
  end
  w = helper1d({v{1}, v{2}, v{3}, v{4}, v{5}, v{6}}, zi, dx, z);
  toc
end



function u = helper1d(f, xi, dx, y)
%1D interpolation helper function.
%
%Interpolate using WENO3 (3 cubic ENO interpolants)
%uses 6 points, y should be between x[3] and x[4] (was 2,3 in C)
%Assumes const dx
%
%There is a maple worksheet that was used to construct this.

  WENOEPS = 1e-6;  % the WENO parameter to prevent div-by-zero

  i = 3;   % was 2 in C code, range below is [i-2,i+3] = [1,6]

  forceChoice = -10;

  % ideal, smooth weights
  C1 = (y-xi-2*dx) .* (-3*dx+y-xi) / (20*dx^2);
  C2 = -(-3*dx+y-xi) .* (y-xi+2*dx) / (10*dx^2);
  C3 = (y-xi+2*dx) .* (y-xi+dx) / (20*dx^2);

  if (forceChoice == -10)
    % compute smoothness indicators
    IS0 = ( -3579.0*f{i+1}.*f{i}  +  2634.0*f{i+1}.*f{i-1}  -  683.0*f{i+1}.*f{i-2}  ...
            -6927.0*f{i}.*f{i-1}  +  1854.0*f{i}.*f{i-2}  -  1659.0*f{i-1}.*f{i-2}  ...
            + 814.0*f{i+1}.*f{i+1}  +  4326.0*f{i}.*f{i}  +  2976.0*f{i-1}.*f{i-1}  ...
            + 244.0*f{i-2}.*f{i-2} ) / 180.0;

    IS1 = ( -3777.0*f{i+1}.*f{i}  +  1074.0*f{i+1}.*f{i-1}  -  1269.0*f{i}.*f{i-1}  ...
            +  1986.0*f{i+1}.*f{i+1}  +  1986.0*f{i}.*f{i}  +  244.0*f{i-1}.*f{i-1}  ...
            +  244.0*f{i+2}.*f{i+2}  -  1269.0*f{i+2}.*f{i+1}  +  1074.0*f{i+2}.*f{i}  ...
            -  293.0*f{i+2}.*f{i-1} ) / 180.0;

    IS2 = ( -3579.0*f{i+1}.*f{i}  +  4326.0*f{i+1}.*f{i+1}  +  814.0*f{i}.*f{i}  ...
            +  2976.0*f{i+2}.*f{i+2}  +  244.0*f{i+3}.*f{i+3}  -  683.0*f{i+3}.*f{i}  ...
            -  6927.0*f{i+2}.*f{i+1}  +  2634.0*f{i+2}.*f{i}  -  1659.0*f{i+3}.*f{i+2}  ...
            +  1854.0*f{i+3}.*f{i+1} ) / 180.0;

    %fprintf('IS0=%g, IS1=%g, IS2=%g\n', IS0, IS1, IS2);

    alpha1 = C1 ./ (WENOEPS + IS0).^2;
    alpha2 = C2 ./ (WENOEPS + IS1).^2;
    alpha3 = C3 ./ (WENOEPS + IS2).^2;
    temp = (alpha1 + alpha2 + alpha3);
    w1 = alpha1 ./ temp;
    w2 = alpha2 ./ temp;
    w3 = alpha3 ./ temp;
  elseif (forceChoice == 10 )
    % forcing ideal weight choice
    w1 = C1;
    w2 = C2;
    w3 = C3;
  elseif (forceChoice == 1)
    % forcing binary right weight choice
    w1 = 0;
    w2 = 0;
    w3 = 1;
  elseif (forceChoice == 0)
    % forcing binary center weight choice
    w1 = 0;
    w2 = 1;
    w3 = 0;
  elseif (forceChoice == -1)
    % forcing binary left weight choice
    w1 = 1;
    w2 = 0;
    w3 = 0;
  end
  % uses f{i-2},...,f{i+1}
  p1 = f{i-2}   -   ( -f{i-1} + f{i-2} ) .* (y-xi+2.0*dx) / dx  ...
       +  ( f{i} - 2.0*f{i-1} + f{i-2} ) .* (y-xi+2.0*dx).*(y-xi+dx) / (2*dx^2)  ...
       -  ( -f{i+1} + 3.0*f{i} - 3.0*f{i-1} + f{i-2} ) .* (y-xi+2.0*dx).*(y-xi+dx).*(y-xi) / (6.0*dx^3);

  % uses f{i-1},...,f{i+2}
  p2 = f{i-1}   +   ( f{i} - f{i-1}) .* (y-xi+dx) / dx  ...
       +  ( f{i+1} - 2.0*f{i} + f{i-1} ) .* (y-xi+dx).*(y-xi) / (2.0*dx^2)  ...
       +  ( f{i+2} - 3.0*f{i+1} + 3.0*f{i} - f{i-1} )  .* (y-xi+dx).*(y-xi).*(y-xi-dx) / (6.0*dx^3);

  % uses f{i},...,f{i+3}
  p3 = f{i}   -   ( -f{i+1} + f{i} ) .* (y-xi) / dx  ...
       +  ( f{i+2} - 2.0*f{i+1} + f{i} ) .* (y-xi).*(y-xi-dx) / (2.0*dx^2)  ...
       -  ( -f{i+3} + 3.0*f{i+2} - 3.0*f{i+1} + f{i} ) .* (y-xi).*(y-xi-dx).*(y-xi-2.0*dx) / (6.0*dx^3);

  %fprintf('w1=%g, p1=%.16g,  w2=%g, p2=%.16g,  w3=%g, p3=%g\n', w1,p1, w2,p2, w3,p3);

  u = w1.*p1 + w2.*p2 + w3.*p3;
end




function Cache = weno6_interp3d_makecache(cp, f, xyz)
%WENO6_INTERP3D_MAKE_CACHE  pre-computation for WENO interpolation 3D

  x1d = cp.x1d;
  y1d = cp.y1d;
  z1d = cp.z1d;
  Nx = length(x1d);
  Ny = length(y1d);
  Nz = length(z1d);

  dx = x1d(2) - x1d(1);  % assumed constant and same in x,y,z

  relpt = cp.x1d(1);  % TODO

  x = xyz(:,1);
  y = xyz(:,2);
  z = xyz(:,3);

  tic
  % determine the basepoint
  [ijk,X] = findGridInterpBasePt(xyz, 5, relpt, dx);
  xi = X(:,1) + 2*dx;
  yi = X(:,2) + 2*dx;
  zi = X(:,3) + 2*dx;
  ijk = ijk + 2;

  I = sub2ind([Ny Nx Nz], ijk(:,2), ijk(:,1), ijk(:,3));

  B = findInBand(I, cp.band, Nx*Ny*Nz);

  [E W N S U D] = neighbourMatrices(cp, cp.band, cp.band);
  toc

  tic
  C = {};
  I = speye(size(W));
  T1 = {W*W, W, I, E, E*E, E*E*E};
  T2 = {S*S, S, I, N, N*N, N*N*N};
  T3 = {D*D, D, I, U, U*U, U*U*U};
  for i=1:6
    for j=1:6
      for k=1:6
        C{i,j,k} = B*T1{i}*T2{j}*T3{k};
      end
    end
  end
  toc

  Cache.C = C;
  Cache.xi = xi;
  Cache.yi = yi;
  Cache.zi = zi;
  Cache.x = x;
  Cache.y = y;
  Cache.z = z;
  Cache.dx = dx;
  Cache.dim = 3;
end




function Cache = weno6_interp2d_makecache(cp, f, xy)
%WENO6_INTERP2D_MAKE_CACHE  pre-computation for WENO interpolation 2D

  x1d = cp.x1d;
  y1d = cp.y1d;
  Nx = length(x1d);
  Ny = length(y1d);

  dx = x1d(2) - x1d(1);  % assumed constant and same in x,y

  relpt = cp.x1d(1);  % TODO

  x = xy(:,1);
  y = xy(:,2);

  tic
  % determine the basepoint
  [ijk,X] = findGridInterpBasePt(xy, 5, relpt, dx);
  xi = X(:,1) + 2*dx;
  yi = X(:,2) + 2*dx;
  ijk = ijk + 2;
  I = sub2ind([Ny Nx], ijk(:,2), ijk(:,1));

  B = findInBand(I, cp.band, Nx*Ny*Nz);

  [E W N S] = neighbourMatrices(cp, cp.band, cp.band);
  toc

  tic
  C = {};
  I = speye(size(W));
  T1 = {W*W, W, I, E, E*E, E*E*E};
  T2 = {S*S, S, I, N, N*N, N*N*N};
  for i=1:6
    for j=1:6
      C{i,j} = B*T1{i}*T2{j};
    end
  end
  toc

  Cache.C = C;
  Cache.xi = xi;
  Cache.yi = yi;
  Cache.x = x;
  Cache.y = y;
  Cache.dx = dx;
  Cache.dim = 2;
end

