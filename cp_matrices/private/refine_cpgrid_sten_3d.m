function gg = refine_cpgrid_sten_3d(g, p)
%REFINE_CPGRID_STEN_3D CP grid refinement in 3D using interp stencil
%
%   Private: you probably want refine_cpgrid()
%
%   gg = refine_cpgrid_sten_3d(cpgrid_old, p)
%     uses (stencils for) degree p interpolation.  p=3 is the
%     default if omitted.

  % see refine_cpgrid_sten_2d for comments

  if nargin < 2
    p = 3;
  end
  assert(g.dim == 3);

  cpfun = g.cpfun;
  dim = g.dim;

  % new grid
  dx = g.dx / 2;
  relpt = [g.x1d(1) g.y1d(1) g.z1d(1)];
  x1d  = ( g.x1d(1):dx:g.x1d(end) )';
  y1d  = ( g.y1d(1):dx:g.y1d(end) )';
  z1d  = ( g.z1d(1):dx:g.z1d(end) )';
  nx = length(x1d);
  ny = length(y1d);
  nz = length(z1d);

  % Main idea
  [Ei,Ej,Es] = interp3_matrix(x1d,y1d,z1d, ...
                              g.cpx,g.cpy,g.cpz, p);
  band = unique(Ej);

  [J,I,K] = ind2sub([ny nx nz], band);

  x = relpt(1) + (I-1)*dx;
  y = relpt(2) + (J-1)*dx;
  z = relpt(3) + (K-1)*dx;

  [cpx,cpy,cpz,dist] = cpfun(xg,yg,zg);


  % Safety check
  [tilde,Ej,tilde] = interp2_matrix(x1d, y1d, cpx, cpy, p);
  band2 = unique(Ej);

  if length(band) == length(band2) && all(band == band2)
    % gave same thing, no need to recompute cp's.
  else
    band = band2;
    [J,I,K] = ind2sub([ny nx nz], band);
    x = relpt(1) + (I-1)*dx;
    y = relpt(2) + (J-1)*dx;
    z = relpt(3) + (K-1)*dx;
    [cpx,cpy,cpz,dist] = cpfun(x,y,z);
  end

  gg.dim = dim;
  gg.dx = dx;
  gg.x1d = x1d;
  gg.y1d = y1d;
  gg.z1d = z1d;
  gg.cpfun = cpfun;
  gg.band = band;
  gg.cpx = cpx;
  gg.cpy = cpy;
  gg.cpz = cpz;
  gg.x = x;
  gg.y = y;
  gg.z = z;
  gg.dist = dist;

