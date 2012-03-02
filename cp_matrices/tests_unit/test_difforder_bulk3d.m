function [pass, str] = test_difforder_bulk3d()
  str = '3D bulk diff matrices: order tests';

  dx = 1/8*pi;
  dy = 1/9*pi;
  dz = 1/10*pi;
  errs1 = helper1(dx,dy,dz);

  dx = dx/2;  dy = dy/2;  dz = dz/2;
  errs2 = helper1(dx,dy,dz);

  orders = errs1 ./ errs2;

  design_ord = [2,2,2, 2,2,2, 1,1,1, 1,1,1, 2,2,2];
  fuzz = 0.7;

  pass = orders > (fuzz * 2.^design_ord);


function err = helper1(dx,dy,dz)
  x1d = -pi:dx:(pi-dx)';
  y1d = -pi:dy:(pi-dy)';
  z1d = -pi:dy:(pi-dz)';
  [x,y,z] = meshgrid(x1d, y1d, z1d);

  [Dxx,Dyy,Dzz, Dxc,Dyc,Dzc, Dxb,Dyb,Dzb, Dxf,Dyf,Dzf, ...
   Dxyc,Dxzc,Dyzc] = bulk3d_matrices(x1d,y1d,z1d);

  xg = x(:);
  yg = y(:);
  zg = z(:);

  u =       sin(2*xg) .*   sin(3*yg) .*   sin(4*zg);
  ux_ex = 2*cos(2*xg) .*   sin(3*yg) .*   sin(4*zg);
  uy_ex =   sin(2*xg) .* (3*cos(3*yg)) .*   sin(4*zg);
  uz_ex =   sin(2*xg) .*   sin(3*yg) .* (4*cos(4*zg));
  uxx_ex = -4*sin(2*xg) .* sin(3*yg) .*   sin(4*zg);
  uyy_ex =    sin(2*xg) .* (-9*sin(3*yg)) .* sin(4*zg);
  uzz_ex =    sin(2*xg) .*    sin(3*yg) .* (-16*sin(4*zg));
  uxy_ex =  2*cos(2*xg) .* (3*cos(3*yg)) .*   sin(4*zg);
  uxz_ex =  2*cos(2*xg) .* sin(3*yg) .*  (4*cos(4*zg));
  uyz_ex =    sin(2*xg) .* (3*cos(3*yg)) .* (4*cos(4*zg));

  uxx = Dxx*u;
  uyy = Dyy*u;
  uzz = Dzz*u;

  ux1 = Dxc*u;
  ux2 = Dxb*u;
  ux3 = Dxf*u;

  uy1 = Dyc*u;
  uy2 = Dyb*u;
  uy3 = Dyf*u;

  uz1 = Dzc*u;
  uz2 = Dzb*u;
  uz3 = Dzf*u;

  uxy = Dxyc*u;
  uxz = Dxzc*u;
  uyz = Dyzc*u;

  c = 0;
  c=c+1; err(c) = max(abs( uxx - uxx_ex ));
  c=c+1; err(c) = max(abs( uyy - uyy_ex ));
  c=c+1; err(c) = max(abs( uzz - uzz_ex ));
  c=c+1; err(c) = max(abs( ux1 - ux_ex ));
  c=c+1; err(c) = max(abs( uy1 - uy_ex ));
  c=c+1; err(c) = max(abs( uz1 - uz_ex ));
  c=c+1; err(c) = max(abs( ux2 - ux_ex ));
  c=c+1; err(c) = max(abs( uy2 - uy_ex ));
  c=c+1; err(c) = max(abs( uz2 - uz_ex ));
  c=c+1; err(c) = max(abs( ux3 - ux_ex ));
  c=c+1; err(c) = max(abs( uy3 - uy_ex ));
  c=c+1; err(c) = max(abs( uz3 - uz_ex ));
  c=c+1; err(c) = max(abs( uxy - uxy_ex ));
  c=c+1; err(c) = max(abs( uxz - uxz_ex ));
  c=c+1; err(c) = max(abs( uyz - uyz_ex ));
