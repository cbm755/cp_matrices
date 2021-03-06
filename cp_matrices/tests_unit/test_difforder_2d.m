function [pass, str] = test_difforder_2d()
  str = 'Order tests of 2D derivatives Dx2c, Dx1b, Dx1f, Dxx2c, etc';

ds = 0.05;
[maxerrx1,maxerry1,maxerrmix1] = helper1(ds);

ds = ds / 2;
[maxerrx2,maxerry2,maxerrmix2] = helper1(ds);

ordersx = maxerrx1 ./ maxerrx2;
ordersy = maxerry1 ./ maxerry2;
ordersmix = maxerrmix1 ./ maxerrmix2;

fuzz = 0.95;

% design orders at 2 1 1 2.
passx = ordersx > (fuzz * 2.^[2 1 1 2]);
passy = ordersy > (fuzz * 2.^[2 1 1 2]);
passm = ordersmix > (fuzz * 2.^2);

pass = [passx passy passm];


function [maxerrx,maxerry,maxerrmix] = helper1(ds)
pad = 4;

R=1;  % Radius
x1d=(-R-pad*ds):ds:(R+pad*ds)';
y1d=(-R-pad*ds):ds:(R+pad*ds)';
nx=length(x1d);
ny=length(y1d);
[x2d,y2d]=meshgrid(x1d,y1d);

[cpx,cpy,dist] = cpCircle(x2d,y2d,R);

dim = 2;
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

xg = x2d(band);
yg = y2d(band);
cpxg = cpx(band);
cpyg = cpy(band);

[Dx1f,Dx1b, Dy1f,Dy1b] = firstderiv_upw1_2d_matrices(x1d,y1d, band,band);
[Dx2c, Dy2c] = firstderiv_cen2_2d_matrices(x1d,y1d, band,band);
[Dxx2c, Dyy2c] = secondderiv_cen2_2d_matrices(x1d,y1d, band,band);
Dxy2c = secondderiv_mixcen2_2d_matrix(x1d,y1d, band,band);

u = sin(2*xg) .* cos(3*yg);
ux_ex = 2*cos(2*xg) .* cos(3*yg);
uy_ex = sin(2*xg) .* (-3*sin(3*yg));
uxx_ex = -4*sin(2*xg) .* cos(3*yg);
uyy_ex = sin(2*xg) .* (-9*cos(3*yg));
uxy_ex = 2*cos(2*xg) .* (-3*sin(3*yg));

ux1 = Dx2c*u;
ux2 = Dx1f*u;
ux3 = Dx1b*u;
uxx = Dxx2c*u;
uy1 = Dy2c*u;
uy2 = Dy1f*u;
uy3 = Dy1b*u;
uyy = Dyy2c*u;
uxy = Dxy2c*u;


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
maxerrmix  = max(abs( uxy(inband)-uxy_ex(inband) ));
if (1==0)
  figure(1); clf; hold on;
  plot(errx, 'r--.')
  plot(ux_ex(inband), 'gx');
  plot(ux1(inband), 'ms');
  
  figure(2); clf; hold on;
  plot(erry, 'r.')
  plot(uy_ex(inband), 'gx');
  plot(uy1(inband), 'ms');
end


% TODO: use polynomials where we know exactly what should happen
