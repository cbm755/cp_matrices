function cpGrid = refine_cpgridnd(cpGrid0, bw)
%REFINE_GRID_BW_ND   Make a finer CP representation based on bandwidth
%
%   Private: you probably want refine_grid_bw()
%
%   cpGrid2 = refine_gridnd(cpGrid, bw)
%

  if ~iscell(cpGrid0.x1d)
    error('require cell arrays for n-D support');
  end

  X1d0 = cpGrid0.x1d;
  dx0 = cpGrid0.dx;
  dist0 = cpGrid0.dist;
  band0 = cpGrid0.band;
  cpf = cpGrid0.cpfun;
  dim = cpGrid0.dim;
  assert(dim == length(X1d0));

  relpt = zeros(1, dim);
  for d = 1:dim
    relpt(d) = X1d0{d}(1);
  end

  Ns0 = zeros(1, dim);
  ddx0 = zeros(1, dim);
  for d=1:dim
    Ns0(d) = length(X1d0{d});
    ddx0(d) = X1d0{d}(2) - X1d0{d}(1);
    if ~assertAlmostEqual(dx0, ddx0(d), 100*eps)
      error('this routine requires the same grid spacing in each direction');
    end
  end

  dx = dx0/2;

  % new 1d grids (the theoretical meshgrid)
  Ns = zeros(1, dim);
  for d = 1:dim
    X1d{d} = (X1d0{d}(1):dx:X1d0{d}(end))';
    Ns(d) = length(X1d{d});
  end

  % The (bw + sqrt(dim)) is to account for using "one-sided"
  % directions below.  This is minimal in the sense that using
  % 0.95*sqrt(dim) changes the results (measured after pruning to
  % the bw).  The other option is to use two-sided differences
  % which then requires a call to "unique" (see below).  This
  % second approach is not yet tested in the n-D case.
  if (1==1)
    % optionally narrow the band to the new bw, could add a safety factor
    I = (abs(dist0) <= (bw + sqrt(dim))*dx);
    band0 = band0(I);
  end

  % neat trick to assign all the outputs to a cell array
  iic = cell(1,dim);
  [iic{:}] = ind2sub(Ns0, band0);


  %% now find those points in the finer grid
  for d=1:dim
    ii0{d} = 2*(iic{d}-1) + 1;
  end
  npts0 = length(ii0{1});



  %% add the neighbours of the finer grid set

  % here we use all 27 neighbours
  dirs = {};

  %ex = cell(1,dim);
  twos = 2*ones(1,d);
  for s=1:prod(twos);
    [ex{1:d}] = ind2sub(twos, s);
    dirs{s} = cell2mat(ex) - 1;
  end
  ndirs = length(dirs);

  %ii = zeros(npts0, dim, ndirs);
  for d = 1:dim
    ii{d} = zeros(npts0, ndirs);
    for n = 1:ndirs
      ii{d}(:, n) = ii0{d} + dirs{n}(d);
    end
  end

  % straighten each
  for d = 1:dim
    ii{d} = ii{d}(:);
  end

  % TODO: if you use non-distinct directions (-1,0,1), must call unique
  %ijk = unique(ijk, 'rows');

  % new band
  band = sub2ind(Ns, ii{:});

  for d = 1:dim
    xg{d} = relpt(d) + (ii{d}-1)*dx;
  end

  [cpx, dist] = cpf(xg);

  % TODO
  %if isempty(bdy0)
  %  [cpx, cpy, cpz, dist] = cpf(xg, yg, zg);
  %  bdy = [];
  %else
  %  [cpx, cpy, cpz, dist, bdy] = cpf(xg, yg, zg);
  %end

  %% our grid should be an overesimate
  I = find(abs(dist) <= bw*dx);
  for d = 1:dim
    xg{d} = xg{d}(I);
    cpx{d} = cpx{d}(I);
  end
  dist = dist(I);
  band = band(I);
  %if ~isempty(bdy0)
  %  bdyg = bdy(I);
  %else
  %  bdyg = [];
  %end

  cpGrid.dim = dim;
  cpGrid.x1d = X1d;
  cpGrid.band = band;
  cpGrid.x = xg;
  cpGrid.cpx = cpx;
  cpGrid.dist = dist;
  cpGrid.cpfun = cpf;
  cpGrid.dx = dx;

