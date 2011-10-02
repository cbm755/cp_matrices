function [IJK,DIST,CP,XYZ] = tri2cp_old(Faces, Vertices, dx, relpt, p, fd_stenrad)
%TRI2CP_OLD  Convert a triangulation to a banded closest point representation
%   Use TRI2CP instead.

%addpath('../../surfaces/tri_pig');
%addpath('../../surfaces/readply');

if (nargin < 4)
  relpt = 0;
end
if (nargin < 5)
  p = 3;  % interpolation degree
end
if (nargin < 6)
  fd_stenrad = 1;  % Finite difference stencil radius
end

dim = 3;

% The formula for bw is found in [Ruuth & Merriman 2008] and the
% 1.00001 is a safety factor.
bw = 1.00001*sqrt((dim-1)*((p+1)/2)^2 + ((fd_stenrad+(p+1)/2)^2));
bw = bw * dx;
%dx = 0.05;


DEBUG_LEVEL = 10;

if (1==0)
  %% load ply file
  %PlyFile = 'pig_loop2.ply';
  PlyFile = 'annies_pig.ply';

  [Faces, Vertices] = plyread(PlyFile, 'tri');
  xp = Vertices(:,1);
  yp = Vertices(:,2);
  zp = Vertices(:,3);


end

%% estimate number of gridpoints (overestimate ok)

%disp('TODO: clean up the hash size estimate');

V = Vertices;

%min1 = min(min(V));
%max1 = max(max(V));
%RAD = max(abs(min1),abs(max1));

%mid = mean(V);
%RAD = max((V(:,1)-mid(1)).^2 + (V(:,2)-mid(2)).^2 + (V(:,3)-mid(3)).^2);
%RAD = sqrt(RAD);

%I like this one
min1 = min(V);
max1 = max(V);
RAD = norm(max1 - min1) / 2;

% We estimate how many grid points we might need based on a sphere,
% then multiply by this much.  Maybe increase if your surface is
% really convoluted!
EXPECTEDGRIDSZ_SAFETY_FACTOR = 6;

% expectedgridsz = ceil((BANDWIDTH*DX)*(2*CONST_PI/DX));
expectedgridsz = EXPECTEDGRIDSZ_SAFETY_FACTOR * ...
    ceil((bw/dx)*(4*pi*RAD^2)/(dx^2));
hashtablesz = ceil(1.25*expectedgridsz);

mex -O CFLAGS='\$CFLAGS -Wall' tri2cp_helper_old.c

tic
[IJK,DD,CP,XYZ] = tri2cp_helper_old(dx, [relpt relpt relpt], bw, ...
                            Faces, Vertices, ...
                            [expectedgridsz, hashtablesz], ...
                            DEBUG_LEVEL);
toc
% NOTE: DD is squared distance
DIST = sqrt(DD);

num_grid_points = length(DD)

% a test
assert( max(max(abs(relpt + (IJK-1)*dx - XYZ))) == 0)
