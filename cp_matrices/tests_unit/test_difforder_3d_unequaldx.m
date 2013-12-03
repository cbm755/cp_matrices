function [pass, str] = test_difforder_3d_unequaldx()
  str = 'Order tests of 3D derivatives and Laplacian, dx != dy != dz';

  ds = 0.25;
  [errx1,erry1,errz1,elap1] = helper1(ds);

  ds = ds / 2;
  [errx2,erry2,errz2,elap2] = helper1(ds);

  ordersx = errx1 ./ errx2;
  ordersy = erry1 ./ erry2;
  ordersz = errz1 ./ errz2;
  ordersL = elap1 ./ elap2;

  fuzz = 0.9;

  % 2 1 1 2 are the design orders
  passx = ordersx > (fuzz * 2.^[2 1 1 2]);
  passy = ordersy > (fuzz * 2.^[2 1 1 2]);
  passz = ordersz > (fuzz * 2.^[2 1 1 2]);
  passL = ordersL > (fuzz * 2.^2);

  pass = [passx passy passz passL];
end


function [maxerrx,maxerry,maxerrz,maxLerr] = helper1(ds)
pad = 5;

% we don't care about a surface here, really just want a banded
% grid to play in
R=1;  % Radius
x1d=(-R-pad*ds):ds:(R+pad*ds)';
y1d=(-R-pad*2*ds):(1.5*ds):(R+pad*2*ds)';
z1d=(-R-pad*ds):(0.9*ds):(R+pad*ds)';
%nx=length(x1d);
%ny=length(y1d);
[x3d,y3d,z3d]=meshgrid(x1d,y1d,z1d);

[cpx,cpy,cpz,dist] = cpSphere(x3d,y3d,z3d,R);

dim = 3;
p = 3;  % degree interp
bw = 1.001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
if (pad < bw)
  warning('maybe insufficient padding');
end
band = find(dist <= bw*ds);
outband = find(dist > bw*ds);
% TODO:testing
inband = find(dist(band) <= (bw/2)*ds);
%inband = find(dist <= (bw/2)*ds);

xg = x3d(band);
yg = y3d(band);
zg = z3d(band);
cpxg = cpx(band);
cpyg = cpy(band);
cpzg = cpz(band);

[Dx1f,Dx1b, Dy1f,Dy1b, Dz1f,Dz1b] = firstderiv_upw1_3d_matrices(x1d,y1d,z1d, band,band);
[Dx2c, Dy2c, Dz2c] = firstderiv_cen2_3d_matrices(x1d,y1d,z1d, band,band);
[Dxx2c, Dyy2c, Dzz2c] = secondderiv_cen2_3d_matrices(x1d,y1d,z1d, band,band);
L = laplacian_3d_matrix(x1d, y1d, z1d, 2, band, band);

u = sin(2*xg) .* cos(3*yg) .* sin(zg);
ux_ex  = 2*cos(2*xg) .* cos(3*yg) .* sin(zg);
uy_ex  = sin(2*xg) .* (-3*sin(3*yg)) .* sin(zg);
uxx_ex = -4*sin(2*xg) .* cos(3*yg) .* sin(zg);
uyy_ex = sin(2*xg) .* (-9*cos(3*yg)) .* sin(zg);
uz_ex  = sin(2*xg) .* cos(3*yg) .* cos(zg);
uzz_ex = sin(2*xg) .* cos(3*yg) .* (-sin(zg));


ux1 = Dx2c*u;
ux2 = Dx1f*u;
ux3 = Dx1b*u;
uxx = Dxx2c*u;
uy1 = Dy2c*u;
uy2 = Dy1f*u;
uy3 = Dy1b*u;
uyy = Dyy2c*u;
uz1 = Dz2c*u;
uz2 = Dz1f*u;
uz3 = Dz1b*u;
uzz = Dzz2c*u;
Lu = L*u;


errx = ux1(inband)-ux_ex(inband);
erry = uy1(inband)-uy_ex(inband);

maxerrx(1) = max(abs( ux1(inband)-ux_ex(inband) ));
maxerrx(2) = max(abs( ux2(inband)-ux_ex(inband) ));
maxerrx(3) = max(abs( ux3(inband)-ux_ex(inband) ));
maxerrx(4) = max(abs( uxx(inband)-uxx_ex(inband) ));
maxerry(1) = max(abs( uy1(inband)-uy_ex(inband) ));
maxerry(2) = max(abs( uy2(inband)-uy_ex(inband) ));
maxerry(3) = max(abs( uy3(inband)-uy_ex(inband) ));
maxerry(4) = max(abs( uyy(inband)-uyy_ex(inband) ));
maxerrz(1) = max(abs( uz1(inband)-uz_ex(inband) ));
maxerrz(2) = max(abs( uz2(inband)-uz_ex(inband) ));
maxerrz(3) = max(abs( uz3(inband)-uz_ex(inband) ));
maxerrz(4) = max(abs( uzz(inband)-uzz_ex(inband) ));
maxLerr = max(abs( Lu(inband) - uxx_ex(inband) ...
                   - uyy_ex(inband) - uzz_ex(inband) ));

end

