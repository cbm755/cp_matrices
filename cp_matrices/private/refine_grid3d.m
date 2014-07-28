function [band,xg,yg,zg,cpxg,cpyg,cpzg,distg,bdyg,dx,x1d,y1d,z1d] = refine_grid3d(cpf,dx0,x1d0,y1d0,z1d0,bw0,band0,dist0,bdy0,use_ndgrid)
%REFINE_GRID3D   Make a finer CP representation
%   Call parent wrapper function "refine_grid_bw" instead.

  dim = 3;

  relpt = [x1d0(1) y1d0(1) z1d0(1)];

  dx = dx0/2;

  Nx0 = length(x1d0);
  Ny0 = length(y1d0);
  Nz0 = length(z1d0);
  % new 1d grids (the theoretical meshgrid)
  x1d = (x1d0(1):dx:x1d0(end))';
  y1d = (y1d0(1):dx:y1d0(end))';
  z1d = (z1d0(1):dx:z1d0(end))';
  Nx = length(x1d);
  Ny = length(y1d);
  Nz = length(z1d);


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
    [iic, jjc, kkc] = ind2sub([Nx0 Ny0 Nz0], band0);
  else
    % meshgrid ordering
    [jjc, iic, kkc] = ind2sub([Ny0 Nx0 Nz0], band0);
  end


  %% now find those points in the finer grid
  ii0 = 2*(iic-1) + 1;
  jj0 = 2*(jjc-1) + 1;
  kk0 = 2*(kkc-1) + 1;
  n = length(ii0);


  %% add the neighbours of the finer grid set

  % here we use all 27 neighbours
  dirs = {};
  d = 0;
  for i=-1:1
    for j=-1:1
      for k=-1:1
        d = d + 1;
        dirs{d} = [i j k];
      end
    end
  end

  ndirs = length(dirs);
  ii = zeros(ndirs*n, 1);
  jj = zeros(ndirs*n, 1);
  kk = zeros(ndirs*n, 1);
  for d=1:ndirs
    ii(1+(d-1)*n:d*n) = ii0 + dirs{d}(1);
    jj(1+(d-1)*n:d*n) = jj0 + dirs{d}(2);
    kk(1+(d-1)*n:d*n) = kk0 + dirs{d}(3);
  end
  ijk = [ii jj kk];
  ijk = unique(ijk, 'rows');
  ii = ijk(:,1);
  jj = ijk(:,2);
  kk = ijk(:,3);

  if use_ndgrid
    warning('not tested');
    band = sub2ind([Nx Ny Nz], ii, jj, kk);
  else
    % meshgrid ordering
    band = sub2ind([Ny Nx Nz], jj, ii, kk);
  end

  xg = relpt(1) + (ii-1)*dx;
  yg = relpt(2) + (jj-1)*dx;
  zg = relpt(3) + (kk-1)*dx;

  if isempty(bdy0)
    [cpx, cpy, cpz, dist] = cpf(xg, yg, zg);
    bdy = [];
  else
    [cpx, cpy, cpz, dist, bdy] = cpf(xg, yg, zg);
  end

  %% our grid should be an overesimate
  I = find(abs(dist) <= bw*dx);
  xg = xg(I);
  yg = yg(I);
  zg = zg(I);
  cpxg = cpx(I);
  cpyg = cpy(I);
  cpzg = cpz(I);
  distg = dist(I);
  band = band(I);
  if ~isempty(bdy0)
    bdyg = bdy(I);
  else
    bdyg = [];
  end
