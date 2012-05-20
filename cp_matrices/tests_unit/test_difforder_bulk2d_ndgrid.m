function [pass, str] = test_difforder_bulk2d_ndgrid()
  str = '2D bulk diff matrices (ndgrid): order tests';

  dx = 0.25*pi/2;
  dy = 0.2*pi/2;
  errs1 = helper1(dx,dy);

  dx = dx/2;  dy = dy/2;
  errs2 = helper1(dx,dy);

  orders = errs1 ./ errs2;

  design_ord = [2,2,2,2,1,1,1,1,2];
  fuzz = 0.9;

  pass = orders > (fuzz * 2.^design_ord);


function err = helper1(dx,dy)
  x1d = -pi:dx:(pi-dx)';
  y1d = -pi:dy:(pi-dy)';
  [x,y] = ndgrid(x1d, y1d);

  [Dxx,Dyy,Dxc,Dyc,Dxb,Dyb,Dxf,Dyf,Dxyc] = bulk2d_matrices(x1d, y1d, true);

  xg = x(:);
  yg = y(:);

  u = sin(2*xg) .* cos(3*yg);
  ux_ex = 2*cos(2*xg) .* cos(3*yg);
  uy_ex = sin(2*xg) .* (-3*sin(3*yg));
  uxx_ex = -4*sin(2*xg) .* cos(3*yg);
  uyy_ex = sin(2*xg) .* (-9*cos(3*yg));
  uxy_ex = 2*cos(2*xg) .* (-3*sin(3*yg));

  ux1 = Dxc*u;
  ux2 = Dxb*u;
  ux3 = Dxf*u;
  uxx = Dxx*u;
  uy1 = Dyc*u;
  uy2 = Dyb*u;
  uy3 = Dyf*u;
  uyy = Dyy*u;
  uxy = Dxyc*u;

  c = 0;
  c=c+1; err(c) = max(abs( uxx - uxx_ex ));
  c=c+1; err(c) = max(abs( uyy - uyy_ex ));
  c=c+1; err(c) = max(abs( ux1 - ux_ex ));
  c=c+1; err(c) = max(abs( uy1 - uy_ex ));
  c=c+1; err(c) = max(abs( ux2 - ux_ex ));
  c=c+1; err(c) = max(abs( uy2 - uy_ex ));
  c=c+1; err(c) = max(abs( ux3 - ux_ex ));
  c=c+1; err(c) = max(abs( uy3 - uy_ex ));
  c=c+1; err(c) = max(abs( uxy - uxy_ex ));
