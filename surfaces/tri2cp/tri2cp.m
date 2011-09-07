addpath('../../surfaces/tri_pig');
addpath('../../surfaces/readply');

dx = 0.1;

relpt = -2;

% use the formula here
p = 3;  bw = 4.1235*dx  % p = 3, 4.12310562561766

%p = 4;  bw = 4.9753*dx  % p = 4, 4.9749371855331

DEBUG_LEVEL = 10;

%% load ply file
%PlyFile = 'pig_loop2.ply';
PlyFile = 'annies_pig.ply';

[Faces, Vertices] = plyread(PlyFile, 'tri');
xp = Vertices(:,1);
yp = Vertices(:,2);
zp = Vertices(:,3);



%% estimate number of gridpoints (overestimate ok)
disp('TODO: clean up the hash size estimate');
V = Vertices;
%min1 = min(min(V));
%max1 = max(max(V));
%RAD = max(abs(min1),abs(max1));

%mid = mean(V);
%RAD = max((V(:,1)-mid(1)).^2 + (V(:,2)-mid(2)).^2 + (V(:,3)-mid(3)).^2);
%RAD = sqrt(RAD);

% I like this one
min1 = min(V);
max1 = max(V);
RAD = norm(max1 - min1) / 2;

% We estimate how many grid points we might need based on a sphere,
% then multiply by this much.  Maybe increase if your surface is
% really convoluted!
EXPECTEDGRIDSZ_SAFETY_FACTOR = 6;

%  //expectedgridsz = ceil((BANDWIDTH*DX)*(2*CONST_PI/DX));
expectedgridsz = EXPECTEDGRIDSZ_SAFETY_FACTOR * ...
    ceil((bw/dx)*(4*pi*RAD^2)/(dx^2))
hashtablesz = ceil(1.25*expectedgridsz)
%dbg_printf(2, "(over)-estimated number of gridpoints: %ld\n", expectedgridsz);
%  dbg_printf(2, "empty hashtable size: %ld\n", hashtablesz);



%% Compile

%mex helper_tri2cp.c
mex -O CFLAGS='\$CFLAGS -Wall -std\=c99' helper_tri2cp.c

tic
[IJK,DD,CP,XYZ] = helper_tri2cp(dx, [relpt, relpt, relpt], bw, ...
                            Faces, Vertices, ...
                            [expectedgridsz, hashtablesz], ...
                            DEBUG_LEVEL);
toc

[IJK(end,:)  DD(end)  CP(end,:)  XYZ(end,:)]

numgp = length(DD);

%assert(numgp == 3594)
%numgp - 3594

% a test
assert( max(max(abs(relpt + (IJK-1)*dx - XYZ))) == 0)

if (1==0)
figure(1); clf;
%uplot = Eplot*u;
%trisurf(Faces,xp,yp,zp, uplot);
trisurf(Faces,xp,yp,zp, zp.^2);
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
shading interp
camlight left
colorbar
pause(0);
end