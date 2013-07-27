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

  [J, I, K] = ind2sub([ny nx nz], band);

  x = relpt(1) + (I-1)*dx;
  y = relpt(2) + (J-1)*dx;
  z = relpt(3) + (K-1)*dx;

  [cpx, cpy, cpz, dist] = cpfun(x, y, z);

  iter = 1;
  while(1)
    [tilde,Ej,tilde] = interp3_matrix(x1d,y1d,z1d, cpx,cpy,cpz, p);
    band2 = unique(Ej);

    if length(band) == length(band2) && all(band == band2)
      fprintf('  refine iter%d: no change, stopping iteration\n', iter);
      break
    else
      [tilde,ii] = setdiff(band2,band);
      [tilde,i4,iii] = intersect(band,band2);
      [tilde,i5] = setdiff(band,band2);
      if isempty(i5)
        fprintf('  refine iter%d: found %d new points\n', iter, length(ii));
      else
        %warning('new band missing entries in first band: can this happen?');
        %length(i5)
        %figure(1); clf; plot3(x,y,z,'r.','markersize',1); hold on;
        %plot3(x(i5),y(i5),z(i5),'ko')
        %keyboard
        fprintf('  refine iter%d: found %d new points (and lost %d others)\n', ...
                iter, length(ii), length(i5));
      end
      [J, I, K] = ind2sub([ny nx nz], band2);   % only for debugging?
      [J2, I2, K2] = ind2sub([ny nx nz], band2(ii));
      % build the xyz coordinates of the points in the band
      x = relpt(1) + (I-1)*dx;
      y = relpt(2) + (J-1)*dx;
      z = relpt(3) + (K-1)*dx;
      xT = relpt(1) + (I2-1)*dx;
      yT = relpt(2) + (J2-1)*dx;
      zT = relpt(3) + (K2-1)*dx;
      % find the closest point
      [cpxT, cpyT, cpzT, distT] = cpfun(xT, yT, zT);
      dist2 = zeros(size(band2));
      dist2(iii) = dist(i4);
      dist2(ii) = distT;
      cpx2 = zeros(size(band2));
      cpx2(iii) = cpx(i4);
      cpx2(ii) = cpxT;
      cpy2 = zeros(size(band2));
      cpy2(iii) = cpy(i4);
      cpy2(ii) = cpyT;
      cpz2 = zeros(size(band2));
      cpz2(iii) = cpz(i4);
      cpz2(ii) = cpzT;
      if (1==0)  % debugging, check it matches recomputing everything
        [cpx3, cpy3, cpz3, dist3] = cpfun(x, y, z);
        assert(norm(dist3-dist2)==0)
        assert(norm(cpx3-cpx2)==0)
        assert(norm(cpy3-cpy2)==0)
        assert(norm(cpz3-cpz2)==0)
      end
      band = band2;
      cpx = cpx2;
      cpy = cpy2;
      cpz = cpz2;
      dist = dist2;
    end
    iter = iter + 1;
  end % while

  fprintf('  refine: found %d points in %d iterations\n', length(band), iter);
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

