function [pass, str] = test_interp_component()
  str = 'interp returns matrix or component, make sure they match (3D)';

  pass = [];

  dx = 0.25;
  x1d = (-2.0:dx:2.0)';
  y1d = x1d;  z1d = x1d;
  nx = length(x1d);
  [xx yy zz] = meshgrid(x1d, y1d, z1d);

  [cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
  cpx = cpx(:); cpy = cpy(:); cpz = cpz(:);

  dim = 3;  p = 3;
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(abs(dist) <= bw*dx);

  % store closest points in the band;
  cpxg = cpx(band); cpyg = cpy(band); cpzg = cpz(band);
  xg = xx(band); yg = yy(band); zg = zz(band);

  % TODO: maybe laplacian should support this too?
  %L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band,band);

  E = interp3_matrix(x1d,y1d,z1d, cpxg, cpyg, cpzg, p, band);
  [Ei,Ej,Es] = interp3_matrix(x1d,y1d,z1d, cpxg, cpyg, cpzg, p, band);

  E2 = sparse(Ei, Ej, Es, length(cpxg), length(cpxg));

  pass = [pass nnz(E-E2)==0];

  % now without banding (full column space)
  E = interp3_matrix(x1d,y1d,z1d, cpxg, cpyg, cpzg, p);
  [Ei,Ej,Es] = interp3_matrix(x1d,y1d,z1d, cpxg, cpyg, cpzg, p);

  E2 = sparse(Ei, Ej, Es, size(E,1), size(E,2));
  pass = [pass nnz(E-E2)==0];

  % now complete matrix
  E = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p);
  [Ei,Ej,Es] = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p);

  E2 = sparse(Ei, Ej, Es, size(E,1), size(E,2));
  pass = [pass nnz(E-E2)==0];

