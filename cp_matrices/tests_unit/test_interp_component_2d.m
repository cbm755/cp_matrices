function [pass, str] = test_interp_component_2d()
  str = 'interp returns matrix or component, make sure they match (2D)';

  pass = [];

  dx = 0.2;
  x1d = (-2.0:dx:2.0)';
  y1d = x1d;
  nx = length(x1d);
  [xx yy] = meshgrid(x1d, y1d);

  [cpx, cpy, dist] = cpCircle(xx,yy);
  cpx = cpx(:); cpy = cpy(:);

  dim = 2;  p = 3;
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(abs(dist) <= bw*dx);

  % store closest points in the band;
  cpxg = cpx(band); cpyg = cpy(band);
  xg = xx(band); yg = yy(band);

  % TODO: maybe laplacian should support this too?
  %L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band,band);

  E = interp2_matrix(x1d, y1d, cpxg, cpyg, p, band);
  [Ei,Ej,Es] = interp2_matrix(x1d, y1d, cpxg, cpyg, p, band);

  E2 = sparse(Ei, Ej, Es, length(cpxg), length(cpxg));

  pass = [pass nnz(E-E2)==0];

  % now without banding (full column space)
  E = interp2_matrix(x1d, y1d, cpxg, cpyg, p);
  [Ei,Ej,Es] = interp2_matrix(x1d, y1d, cpxg, cpyg, p);

  E2 = sparse(Ei, Ej, Es, size(E,1), size(E,2));
  pass = [pass nnz(E-E2)==0];

  % now complete matrix
  E = interp2_matrix(x1d, y1d, cpx, cpy, p);
  [Ei,Ej,Es] = interp2_matrix(x1d, y1d, cpx, cpy, p);

  E2 = sparse(Ei, Ej, Es, size(E,1), size(E,2));
  pass = [pass nnz(E-E2)==0];

