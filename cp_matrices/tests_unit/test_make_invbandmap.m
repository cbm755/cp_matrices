function [pass, str] = test_make_invbandmap()
  str = 'test inverse band map';

c = 0;
pass = [];

%% 2D
x1d = linspace(-2, 2.3, 10);
y1d = linspace(-3, 2.1, 12);

[xx,yy] = meshgrid(x1d, y1d);

xx = xx(:);
yy = yy(:);

band = unique(ceil(length(xx)*rand(20,1)));

M = length(x1d)*length(y1d);
invband = make_invbandmap(M, band);

c = c + 1;
pass(c) = invband(band(1)) == 1;

c = c + 1;
pass(c) = all(invband(band) == (1:length(band))');

% doesn't really test the inverse map, just a sanity check
xg = xx(band);
yg = yy(band);
b = band(end);
c = c + 1; pass(c) = xg(end) == xx(b);
c = c + 1; pass(c) = yg(end) == yy(b);
[j,i] = ind2sub([12,10], b);
c = c + 1; pass(c) = x1d(i) == xx(b);
c = c + 1; pass(c) = y1d(j) == yy(b);



%% 2D using cpgrid
clear cpgrid
cpgrid.dim = 2;
cpgrid.x1d = x1d;
cpgrid.y1d = y1d;
cpgrid.band = band;

invband = make_invbandmap(cpgrid);
c = c + 1;
pass(c) = all(invband(cpgrid.band) == (1:length(cpgrid.band))');


%% 2D using cpgrid with n-D form
clear cpgrid
cpgrid.dim = 2;
cpgrid.x1d = {x1d y1d};
cpgrid.band = band;

invband = make_invbandmap(cpgrid);
c = c + 1;
pass(c) = all(invband(cpgrid.band) == (1:length(cpgrid.band))');


%% 3D
x1d = linspace(-2, 2.3, 10);
y1d = linspace(-3, 2.1, 12);
z1d = linspace(-1.2, 1.2, 2);

[xx,yy,zz] = meshgrid(x1d, y1d, z1d);

xx = xx(:);
yy = yy(:);
zz = zz(:);

band = unique(ceil(length(xx)*rand(20,1)));

M = length(x1d)*length(y1d)*length(z1d);
invband = make_invbandmap(M, band);

c = c + 1;
pass(c) = all(invband(band) == (1:length(band))');



%% 3D using cpgrid
clear cpgrid
cpgrid.dim = 3;
cpgrid.x1d = x1d;
cpgrid.y1d = y1d;
cpgrid.z1d = z1d;
cpgrid.band = band;

invband = make_invbandmap(cpgrid);
c = c + 1;
pass(c) = all(invband(cpgrid.band) == (1:length(cpgrid.band))');



%% 3D using cpgrid, n-D form
clear cpgrid
cpgrid.dim = 3;
cpgrid.x1d = {x1d y1d z1d};
cpgrid.band = band;

invband = make_invbandmap(cpgrid);
c = c + 1;
pass(c) = all(invband(cpgrid.band) == (1:length(cpgrid.band))');


%% 4D test for n-D
x1d = linspace(-2, 2.3, 5);
y1d = linspace(-3, 2.1, 6);
z1d = linspace(-1.2, 1.2, 7);
w1d = linspace(-1, 1, 8);

clear xx;

[xx{1:4}] = ndgrid(x1d, y1d, z1d, w1d);

band = unique(ceil(5*6*7*8*rand(20,1)));

clear cpgrid
cpgrid.dim = 4;
cpgrid.x1d = {x1d y1d z1d w1d};
cpgrid.band = band;

invband = make_invbandmap(cpgrid);

c = c + 1;
pass(c) = all(invband(band) == (1:length(band))');


%% 4D using NN form
M = length(x1d)*length(y1d)*length(z1d)*length(w1d);
invband = make_invbandmap(M, band);
c = c + 1;
pass(c) = all(invband(band) == (1:length(band))');

