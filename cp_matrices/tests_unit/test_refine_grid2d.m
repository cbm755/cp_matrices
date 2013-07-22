function [pass, str] = test_refine_grid2d()
  str = 'refinement test: compare band from fine grid with refined grid';

  xlim = [-2.2 2.2];
  ylim = [-1.6 1.6];
  cpf = @cpEllipse;
  %cpf = @cpCircle;

  % standard banding forumla
  dim = 2;
  p = 3;
  bw = rm_bandwidth(dim, p);

  tic
  [band,x,y,cpx,cpy,dist,dx,x1d,y1d] = make_coarse_and_refine(cpf, xlim, ylim, bw, .2, 2);
  time_refine = toc
  tic
  [band2,dx2] = make_fine(cpf, xlim, ylim, bw, 0.1);
  time_fine = toc

  passdx = dx == dx2;
  passlen = length(band) == length(band2);
  passband = max(abs(sort(band)-sort(band2))) == 0;

  [j,i] = ind2sub([length(y1d) length(x1d)], band);
  passxcoor = all( abs(x1d(1)+dx*(i-1)  -  x) < 5*eps);
  passycoor = all( abs(y1d(1)+dx*(j-1)  -  y) < 5*eps);
  calcdist = sqrt( (x - cpx).^2 + (y - cpy).^2 );
  passdist = all( abs(abs(dist) - calcdist) < 5*eps);
  pass = [passdx passxcoor passycoor passdist passlen passband];
  %keyboard
end


function [band,dx] = make_fine(cpf, xlim, ylim, bw, dx)
  x1d = (xlim(1):dx:xlim(2))';
  y1d = (ylim(1):dx:ylim(2))';

  [xx yy] = meshgrid(x1d, y1d);
  [cpx, cpy, dist] = cpf(xx,yy);

  band = find(abs(dist) <= bw*dx);
end


function [band, xg, yg, cpxg, cpyg, dist, dx, x1d, y1d] = ...
      make_coarse_and_refine(cpf, xlim, ylim, bw, dx, NLevels)

  x1d = (xlim(1):dx:xlim(2))';
  y1d = (ylim(1):dx:ylim(2))';

  [xx yy] = meshgrid(x1d, y1d);
  [cpx, cpy, dist] = cpf(xx,yy);

  band = find(abs(dist) <= bw*dx);

  % store closest points in the band (can discard others)
  cpxg = cpx(band); cpyg = cpy(band);
  xg = xx(band); yg = yy(band);
  distg = dist(band);

  [band, xg,yg, cpxg,cpyg, dist, dx, x1d, y1d] = ...
      refine_grid(NLevels-1, cpf, dx, x1d, y1d, bw, band, distg);
end

