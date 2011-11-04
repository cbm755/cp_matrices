function [pass, str] = test_refine_grid()
  str = 'refinment test: compare band from fine grid with refined grid';

  % Include the cp_matrices folder (edit as appropriate)
  addpath('../../cp_matrices');

  % add functions for finding the closest points
  addpath('../../surfaces');

  xlim = [-2.2 2.2];
  ylim = [-1.6 1.6];
  cpf = @cpEllipse;
  %cpf = @cpCircle;

  % standard banding forumla
  dim = 2;
  p = 3;
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));

  tic
  [dx1,band1] = make_coarse_and_refine(cpf, xlim, ylim, bw, .2, 3);
  time_refine = toc
  tic
  [dx2,band2] = make_fine(cpf, xlim, ylim, bw, 0.05);
  time_fine = toc

  passdx = dx1 == dx2;
  passlen = length(band1) == length(band2);
  passband = max(abs(sort(band1)-sort(band2))) ~= 0;
  pass = [passdx passlen passband];
end


function [dx,band] = make_fine(cpf, xlim, ylim, bw, dx)

  x1d = (xlim(1):dx:xlim(2))';
  y1d = (ylim(1):dx:ylim(2))';

  [xx yy] = meshgrid(x1d, y1d);
  [cpx, cpy, dist] = cpf(xx,yy);

  band = find(abs(dist) <= bw*dx);
end


function [mydx, myband] = make_coarse_and_refine(cpf, xlim, ylim, bw, ...
                                             dx, NLevels)

  x1d = (xlim(1):dx:xlim(2))';
  y1d = (ylim(1):dx:ylim(2))';

  relpt = [x1d(1) y1d(1)];

  nx = length(x1d);
  ny = length(y1d);

  [xx yy] = meshgrid(x1d, y1d);
  [cpx, cpy, dist] = cpf(xx,yy);

  band = find(abs(dist) <= bw*dx);

  % store closest points in the band (can discard others)
  cpxg = cpx(band); cpyg = cpy(band);
  xg = xx(band); yg = yy(band);
  distg = dist(band);


  % refine, first level is what we already computed
  k = 1;
  a_band{k} = band;
  a_xg{k} = xg;
  a_yg{k} = yg;
  a_cpx{k} = cpxg;
  a_cpy{k} = cpyg;
  a_dist{k} = distg;
  a_x1d{k} = x1d;
  a_y1d{k} = y1d;
  a_dx{k} = dx;

  % loop to refine
  for k=2:NLevels
    a_dx{k} = a_dx{k-1}/2;
    [a_band{k}, a_xg{k}, a_yg{k}, a_cpx{k}, a_cpy{k}, a_dist{k}, a_x1d{k}, a_y1d{k}] = ...
        refine_grid(cpf, a_dx{k-1}, a_x1d{k-1}, a_y1d{k-1}, bw, a_band{k-1}, a_dist{k-1});
  end

  myband = a_band{end};
  mydx = a_dx{end};
end