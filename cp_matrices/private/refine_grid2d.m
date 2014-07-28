function [band,xg,yg,cpxg,cpyg,distg,bdyg,dx,x1d,y1d] = refine_grid2d(cpf,dx0,x1d0,y1d0,bw0,band0,dist0,bdy0,use_ndgrid)
%REFINE_GRID2D   Make a finer CP representation
%   Call parent wrapper function "refine_grid_bw" instead.

  dim = 2;

  relpt = [x1d0(1) y1d0(1)];

  dx = dx0/2;

  Nx0 = length(x1d0);
  Ny0 = length(y1d0);
  % new 1d grids (the theoretical meshgrid)
  x1d = (x1d0(1):dx:x1d0(end))';
  y1d = (y1d0(1):dx:y1d0(end))';
  Nx = length(x1d);
  Ny = length(y1d);


  bw = bw0;  % TODO: could recalculate (at least give the user and
              % option to specify a new bandwidth

  if (1==1)
    % optionally narrow the band to the new bw, could add a safety factor
    I = (abs(dist0) <= (bw + 0*sqrt(2))*dx);
    band0 = band0(I);
  end


  if use_ndgrid
    % ndgrid ordering: not tested!
    warning('ndgrid not tested');
    [iic, jjc] = ind2sub([Nx0 Ny0], band0);
  else
    % meshgrid ordering
    [jjc, iic] = ind2sub([Ny0 Nx0], band0);
  end


  %% now find those points in the finer grid
  ii0 = 2*(iic-1) + 1;
  jj0 = 2*(jjc-1) + 1;
  n = length(ii0);


  %% add the neighbours of the finer grid set

  % the 9 neighbours including the point itself
  dirs = {[-1 -1] [0 -1] [1 -1] [-1 0] [0 0] [1 0] [-1 1] [0 1] [1 1] };
  % simplier but would bias to upper-left
  %dirs = {[0 0] [1 0] [0 1] [1 1]};
  ndirs = length(dirs);
  ii = zeros(ndirs*n, 1);
  jj = zeros(ndirs*n, 1);
  for d=1:ndirs
    ii(1+(d-1)*n:d*n) = ii0 + dirs{d}(1);
    jj(1+(d-1)*n:d*n) = jj0 + dirs{d}(2);
  end
  ij = [ii jj];
  ij = unique(ij,'rows');
  ii = ij(:,1);
  jj = ij(:,2);

  if use_ndgrid
    warning('not tested');
    band = sub2ind([Nx Ny], ii, jj);
  else
    % meshgrid ordering
    band = sub2ind([Ny Nx], jj, ii);
  end

  xg = relpt(1) + (ii-1)*dx;
  yg = relpt(2) + (jj-1)*dx;

  if isempty(bdy0)
    [cpx, cpy, dist] = cpf(xg, yg);
    bdy = [];
  else
    [cpx, cpy, dist, bdy] = cpf(xg, yg);
  end

  %% our grid should be an overesimate
  I = find(abs(dist) <= bw*dx);
  xg = xg(I);
  yg = yg(I);
  cpxg = cpx(I);
  cpyg = cpy(I);
  distg = dist(I);
  band = band(I);
  if ~isempty(bdy0)
    bdyg = bdy(I);
  else
    bdyg = [];
  end


  makePlots = 0;
  if (makePlots)
    figure(10); clf;
    figure(11); clf;
    figure(12); clf;
    u0 = sin(xg0);
    plot2d_compdomain(u0,xg0,yg0,dx0,dx0,10)
    %axis([-1.5 1.5 -1.5 1.5]);
    %plot(xg0,yg0,'k.');
    xp = relpt(1) + (ii0-1)*dx;
    yp = relpt(2) + (jj0-1)*dx;
    plot(xp, yp, 'k.');
    z=exp(1i*(0:(2*pi/100):2*pi));
    plot(real(z),imag(z),'k-', 'linewidth', 2);

    u = sin(xg);
    plot2d_compdomain(u,xg,yg,dx,dx,11)
    %axis([-1.5 1.5 -1.5 1.5]);
    plot(xp, yp, 'k.');
    plot(real(z),imag(z),'k-', 'linewidth', 2);

    u = sin(xg);
    plot2d_compdomain(u,xg,yg,dx,dx,12)
    %axis([-1.5 1.5 -1.5 1.5]);
    plot(xp,yp,'k.');
    plot(real(z),imag(z),'k-', 'linewidth', 2);
  end
