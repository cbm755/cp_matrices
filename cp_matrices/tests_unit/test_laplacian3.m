function [pass, str] = test_laplacian3()
  str = '3D laplacian matches component-by-component';

  ds = 0.1;

  pad = 5;
  R=1;  % Radius
  x1d=(-R-pad*ds):ds:(R+pad*ds)';
  y1d=(-R-(pad+1)*ds):ds:(R+(pad+1)*ds)';
  z1d=(-R-(pad+2)*ds):ds:(R+(pad+2)*ds)';
  [x,y,z]=meshgrid(x1d,y1d,z1d);

  [cpx,cpy,cpz,dist] = cpSphere(x,y,z,R);

  dim = 3;
  p = 3;  % degree interp
  order = 2;  % laplacian order
  bw = 1.001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(dist <= bw*ds);

  xg = x(band);
  yg = y(band);
  zg = z(band);
  cpxg = cpx(band);
  cpyg = cpy(band);
  cpzg = cpz(band);

  L = laplacian_3d_matrix(x1d,y1d,z1d, order, band,band);

  [Dxx2c, Dyy2c, Dzz2c] = secondderiv_cen2_3d_matrices(x1d,y1d,z1d, band,band);

  u = sin(2*xg) .* cos(3*yg) .* cos(zg);

  uxx = Dxx2c*u;
  uyy = Dyy2c*u;
  uzz = Dzz2c*u;

  lapu = L*u;

  err = max(abs(uxx+uyy+uzz-lapu));

  pass = assertAlmostEqual(err, 0, 1000*eps);

