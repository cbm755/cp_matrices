function gg = cpgrid_restrict_bw(g, bw)
%CPGRID_RESTRICT_BW
%
%  g2 = cpgrid_restrict_bw(g, bw)

  if iscell(g.x1d)
    %% General ndgrid case
    gg.dim = g.dim;
    gg.dx = g.dx;
    gg.x1d = g.x1d;
    gg.cpfun = g.cpfun;
    band = find(abs(g.dist) <= bw*g.dx);
    gg.band = band;
    for d=1:g.dim
      gg.x{d} = g.x{d}(band);
      gg.cpx{d} = g.cpx{d}(band);
    end
    gg.dist = g.dist(band);
  elseif g.dim == 2
    %% 2D
    gg.dim = g.dim;
    gg.dx = g.dx;
    gg.x1d = g.x1d;
    gg.y1d = g.y1d;
    gg.cpfun = g.cpfun;
    band = find(abs(g.dist) <= bw*g.dx);
    gg.band = band;
    gg.x = g.x(band);
    gg.y = g.y(band);
    gg.cpx = g.cpx(band);
    gg.cpy = g.cpy(band);
    gg.dist = g.dist(band);
  elseif g.dim == 3
    %% 3D
    gg.dim = g.dim;
    gg.dx = g.dx;
    gg.x1d = g.x1d;
    gg.y1d = g.y1d;
    gg.z1d = g.z1d;
    gg.cpfun = g.cpfun;
    band = find(abs(g.dist) <= bw*g.dx);
    gg.band = band;
    gg.x = g.x(band);
    gg.y = g.y(band);
    gg.z = g.z(band);
    gg.cpx = g.cpx(band);
    gg.cpy = g.cpy(band);
    gg.cpz = g.cpz(band);
    gg.dist = g.dist(band);
  else
    error('this grid object cannot be handled');
  end


