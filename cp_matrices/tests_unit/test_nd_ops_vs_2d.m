function [pass, str] = test_nd_ops_vs_2d()
  str = 'n-D matrices test: compare to output from 2D code';

  pass = [];
  c = 0;

  dx = 0.1;

  % make vectors of x, y, positions of the grid
  x1d = (-2.0:dx:2.0)';
  y1d = (-2.1:dx:2.1)';

  nx = length(x1d);
  ny = length(y1d);

  [xx yy] = ndgrid(x1d, y1d);

  [cpx,cpy,dist] = cpEllipse(xx,yy, 1.2, 0.9, dx*[1/3 1/4]);

  % banding
  dim = 2;
  p = 3;       % max interpolation order
  stenrad = 2; % max stencil radius for finite differences
  bw = rm_bandwidth(dim, p, stenrad);
  band = find(abs(dist) <= bw*dx);

  % store closest points in the band;
  cpx = cpx(band); cpy = cpy(band);
  xg = xx(band); yg = yy(band);


  % interpolation

  E = interp2_matrix(x1d,y1d,  cpx, cpy, 3, band, true);
  En = interpn_matrix({x1d y1d}, {cpx cpy}, 3, band);

  c = c + 1;
  pass(c) = nnz(E - En) == 0;

  E = interp2_matrix(x1d,y1d,  cpx, cpy, 1, band, true);
  En = interpn_matrix({x1d y1d}, {cpx cpy}, 1, band);

  c = c + 1;
  pass(c) = nnz(E - En) == 0;


  % first deriv: upwinding

  [Dxb,Dxf, Dyb,Dyf] = ...
      firstderiv_upw1_2d_matrices(x1d,y1d, band, band, true);

  [Db,Df] = firstderiv_upw1_nd_matrices({x1d,y1d}, band);

  c = c + 1;  pass(c) = nnz(Dxb - Db{1}) == 0;
  c = c + 1;  pass(c) = nnz(Dyb - Db{2}) == 0;
  c = c + 1;  pass(c) = nnz(Dxf - Df{1}) == 0;
  c = c + 1;  pass(c) = nnz(Dyf - Df{2}) == 0;


  % first deriv: centered

  [Dxc,Dyc] = firstderiv_cen2_2d_matrices(x1d,y1d, band, band, true);
  Dc = firstderiv_cen2_nd_matrices({x1d,y1d}, band);

  c = c + 1;  pass(c) = nnz(Dxc - Dc{1}) == 0;
  c = c + 1;  pass(c) = nnz(Dyc - Dc{2}) == 0;


  % second deriv

  [Dxxc,Dyyc] = secondderiv_cen2_2d_matrices(x1d,y1d, band, band, true);
  D2c = secondderiv_cen2_nd_matrices({x1d,y1d}, band);

  c = c + 1;  pass(c) = nnz(Dxxc - D2c{1}) == 0;
  c = c + 1;  pass(c) = nnz(Dyyc - D2c{2}) == 0;


  % Laplacians

  L2 = laplacian_2d_matrix(x1d,y1d, 2, band, band, true);
  L2n = laplacian_nd_matrix({x1d,y1d}, 2, band, band);

  c = c + 1;
  pass(c) = nnz(L2 - (D2c{1} + D2c{2})) == 0;
  c = c + 1;
  pass(c) = nnz(L2 - L2n) == 0;

  L4 = laplacian_2d_matrix(x1d,y1d, 4, band, band, true);
  L4n = laplacian_nd_matrix({x1d,y1d}, 4, band, band);

  c = c + 1;
  pass(c) = nnz(L4 - L4n) == 0;

