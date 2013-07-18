function [pass, str] = test_difforder_mixed_3d()
  str = 'Order tests of 3D mixed partials';

ds = 0.25;
[err1] = helper1(ds);

ds = ds / 2;
[err2] = helper1(ds);

orders = err1 ./ err2;

fuzz = 0.9;

% [ 2 2 2 ] is design order
passvec = orders > (fuzz * 2.^[2 2 2]);

pass = all(passvec);


function [maxerr] = helper1(ds)
pad = 5;

R=1;  % Radius
x1d=(-R-pad*ds):ds:(R+pad*ds)';
y1d=(-R-pad*ds):ds:(R+pad*ds)';
z1d=(-R-pad*ds):ds:(R+pad*ds)';
nx=length(x1d);
ny=length(y1d);
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

[Dxy2c, Dxz2c, Dyz2c] = secondderiv_mixcen2_3d_matrices(x1d,y1d,z1d, band,band);

u = sin(2*xg) .* cos(3*yg) .* sin(zg);
uxy_ex = 2*cos(2*xg) .* (-3*sin(3*yg)) .* sin(zg);
uxz_ex = 2*cos(2*xg) .* cos(3*yg) .* cos(zg);
uyz_ex = sin(2*xg) .* (-3*sin(3*yg)) .* cos(zg);

uxy = Dxy2c*u;
uxz = Dxz2c*u;
uyz = Dyz2c*u;



maxerr(1) = max(abs( uxy(inband)-uxy_ex(inband) ));
maxerr(2) = max(abs( uxz(inband)-uxz_ex(inband) ));
maxerr(3) = max(abs( uyz(inband)-uyz_ex(inband) ));
