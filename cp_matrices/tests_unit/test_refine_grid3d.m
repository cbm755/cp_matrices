function [pass, str] = test_refine_grid3d()
  str = 'refinement test: compare band from fine grid with refined grid';

  xlim = [-2.2 2.2];
  ylim = [-2.2 2.2];
  zlim = [-1.2 2.2];
  cpf = @cpHemisphere;
  % cpf = cpbar

  % standard banding forumla
  dim = 3;
  p = 3;
  bw = rm_bandwidth(dim, p);

  tic
  [band,x,y,z,cpx,cpy,cpz,dist,bdy,dx,x1d,y1d,z1d] = ...
      make_coarse_and_refine(cpf, xlim, ylim, zlim, bw, .2, 3);
  time_refine = toc
  tic
  [band2,dx2] = make_fine(cpf, xlim, ylim, zlim, bw, 0.05);
  time_fine = toc

  passdx = dx == dx2;
  passlen = length(band) == length(band2);
  passband = max(abs(sort(band)-sort(band2))) == 0;

  [j,i,k] = ind2sub([length(y1d) length(x1d) length(z1d)], band);
  passxcoor = all( abs(x1d(1)+dx*(i-1)  -  x) < 5*eps);
  passycoor = all( abs(y1d(1)+dx*(j-1)  -  y) < 5*eps);
  passzcoor = all( abs(z1d(1)+dx*(k-1)  -  z) < 5*eps);
  calcdist = sqrt( (x - cpx).^2 + (y - cpy).^2 + (z - cpz).^2 );
  passdist = all( abs(abs(dist) - calcdist) < 5*eps);
  pass = [passdx passxcoor passycoor passzcoor passdist passlen passband];
end


function [band,dx] = make_fine(cpf, xlim, ylim, zlim, bw, dx)
  x1d = (xlim(1):dx:xlim(2))';
  y1d = (ylim(1):dx:ylim(2))';
  z1d = (zlim(1):dx:zlim(2))';

  [xx yy zz] = meshgrid(x1d, y1d, z1d);
  [cpx, cpy, cpz, dist] = cpf(xx,yy,zz);

  band = find(abs(dist) <= bw*dx);
end


function [band, x,y,z, cpx,cpy,cpz, dist, bdy, dx, x1d,y1d,z1d] = ...
      make_coarse_and_refine(cpf, xlim, ylim, zlim, bw, dx, NLevels)

  x1d = (xlim(1):dx:xlim(2))';
  y1d = (ylim(1):dx:ylim(2))';
  z1d = (zlim(1):dx:zlim(2))';

  [xx yy zz] = meshgrid(x1d, y1d, z1d);
  [cpx, cpy, cpz, dist, bdy] = cpf(xx,yy,zz);

  band = find(abs(dist) <= bw*dx);

  % store closest points in the band (can discard others)
  cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
  x = xx(band); y = yy(band); z = zz(band);
  dist = dist(band);
  bdy = bdy(band);

  [band, x,y,z, cpx,cpy,cpz, dist, bdy, dx, x1d,y1d,z1d] = ...
      refine_grid(NLevels-1, cpf, dx, x1d,y1d,z1d, bw, band, dist, bdy);
end

