function [pass, str] = test_laplacian2()
  str = '2D laplacian matches component-by-component';

  ds = 0.1;

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
order = 2;  % laplacian order
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

L = laplacian_2d_matrix(x1d,y1d, order, band,band);

[Dxx2c, Dyy2c] = secondderiv_cen2_2d_matrices(x1d,y1d, band,band);

u = sin(2*xg) .* cos(3*yg);
% TODO: or random

uxx = Dxx2c*u;

uyy = Dyy2c*u;

lapu = L*u;

err = max(abs(uxx+uyy-lapu));

pass = assertAlmostEqual(err, 0, 1000*eps);
