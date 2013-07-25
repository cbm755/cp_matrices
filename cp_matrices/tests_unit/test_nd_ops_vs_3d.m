function [pass, str] = test_nd_ops_vs_3d()
  str = 'n-D matrices test: compare to output from 3D code';

  pass = [];
  c = 0;

  dx = 0.1;

  % make vectors of x, y, positions of the grid
  x1d = (-2.0:dx:2.0)';
  y1d = (-2.1:dx:2.1)';
  z1d = (-2.2:dx:2.2)';

  nx = length(x1d);
  ny = length(y1d);
  nz = length(z1d);

  [xx yy zz] = ndgrid(x1d, y1d, z1d);

  [cpx,cpy,cpz,dist] = cpSphere(xx,yy,zz, 1, dx*[1/3 1/4 1/5]);


  % banding
  dim = 3;
  p = 3;       % max interpolation order
  stenrad = 2; % max stencil radius for finite differences
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((stenrad+(p+1)/2)^2));
  band = find(abs(dist) <= bw*dx);

  % store closest points in the band;
  cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
  xg = xx(band); yg = yy(band); zg = zz(band);


  % interpolation

  E = interp3_matrix(x1d,y1d,z1d,  cpx, cpy, cpz, 3, band, true);
  En = interpn_matrix({x1d y1d z1d}, [cpx cpy cpz], 3, band);
  c = c + 1;
  pass(c) = nnz(E - En) == 0;

  % same with cell-array
  En2 = interpn_matrix({x1d y1d z1d}, {cpx cpy cpz}, 3, band);
  c = c + 1;
  pass(c) = nnz(En - En2) == 0;

  E = interp3_matrix(x1d,y1d,z1d,  cpx, cpy, cpz, 1, band, true);
  En = interpn_matrix({x1d y1d z1d}, {cpx cpy cpz}, 1, band);

  c = c + 1;
  pass(c) = nnz(E - En) == 0;


  % first deriv: upwinding

  [Dxb,Dxf, Dyb,Dyf, Dzb,Dzf] = ...
      firstderiv_upw1_3d_matrices(x1d,y1d,z1d, band, band, true);

  [Db,Df] = firstderiv_upw1_nd_matrices({x1d,y1d,z1d}, band);

  c = c + 1;  pass(c) = nnz(Dxb - Db{1}) == 0;
  c = c + 1;  pass(c) = nnz(Dyb - Db{2}) == 0;
  c = c + 1;  pass(c) = nnz(Dzb - Db{3}) == 0;
  c = c + 1;  pass(c) = nnz(Dxf - Df{1}) == 0;
  c = c + 1;  pass(c) = nnz(Dyf - Df{2}) == 0;
  c = c + 1;  pass(c) = nnz(Dzf - Df{3}) == 0;


  % first deriv: centered

  [Dxc,Dyc,Dzc] = firstderiv_cen2_3d_matrices(x1d,y1d,z1d, band, band, true);
  Dc = firstderiv_cen2_nd_matrices({x1d,y1d,z1d}, band);

  c = c + 1;  pass(c) = nnz(Dxc - Dc{1}) == 0;
  c = c + 1;  pass(c) = nnz(Dyc - Dc{2}) == 0;
  c = c + 1;  pass(c) = nnz(Dzc - Dc{3}) == 0;


  % second deriv

  [Dxxc,Dyyc,Dzzc] = secondderiv_cen2_3d_matrices(x1d,y1d,z1d, band, band, true);
  D2c = secondderiv_cen2_nd_matrices({x1d,y1d,z1d}, band);

  c = c + 1;  pass(c) = nnz(Dxxc - D2c{1}) == 0;
  c = c + 1;  pass(c) = nnz(Dyyc - D2c{2}) == 0;
  c = c + 1;  pass(c) = nnz(Dzzc - D2c{3}) == 0;


  % Laplacians

  L2 = laplacian_3d_matrix(x1d,y1d,z1d, 2, band, band, true);
  L2n = laplacian_nd_matrix({x1d,y1d,z1d}, 2, band, band);

  c = c + 1;
  pass(c) = nnz(L2 - (D2c{1} + D2c{2} + D2c{3})) == 0;
  c = c + 1;
  pass(c) = nnz(L2 - L2n) == 0;

  L4 = laplacian_3d_matrix(x1d,y1d,z1d, 4, band, band, true);
  L4n = laplacian_nd_matrix({x1d,y1d,z1d}, 4, band, band);

  c = c + 1;
  pass(c) = nnz(L4 - L4n) == 0;

