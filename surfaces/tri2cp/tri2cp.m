function [IJK,DIST,CP,XYZ,CPFACE] = tri2cp(Faces, Vertices, dx, relpt, p, fd_stenrad)
%TRI2CP  Convert a triangulation to a banded closest point representation
%   Given a triangulation as a list of faces and vertices (for
%   example, from a .ply file), TRI2CP will convert them to a
%   closest point representation.
%
%   [IJK,DIST,CP,XYZ] = tri2cp(Faces, Vertices, dx, relpt, p, fd_stenrad)
%
%   'Faces', 'Vertices': a triangulation (e.g., see triplot).  Note
%                        the 'Faces' should be indexed Matlab-style
%                        (from 1 not 0).
%
%   'dx': the grid size.
%
%   'relpt': a 3-vector specifying a point from which the integer grid
%            is defined.  A scalar relpt will be promoted to a
%            3-vector.  (if you are debugging and using the full 3D
%            matrix support, then relpt must be lower left corner of
%            data.)
%
%   'p': degree of interpolation
%
%   'fd_stenrad': the "radius" of the finite difference stencil (this
%                 is 1 for the 2nd-order Laplacian, 2 for the
%                 4th-order Laplacian).
%
%   [IJK,DIST,CP,XYZ] = tri2cp(Faces, Vertices, dx, relpt, bandwidth)
%
%   'bandwidth': Alternatively, specify the bandwidth directly instead
%                of calculating it from 'p' and 'fd_stenrad'.
%
%   *** OUTPUT ***
%
%   'IJK': a 3 column matrix, each row of which gives the indices of a
%          gridpt relative to 'relpt'.  These are Matlab-style (i.e.,
%          i=1 corresponds to relpt).
%
%   'DIST': distance from each gridpt to its closest point.
%
%   'CP': a 3 column matrix, the closest point of each gridpt.
%
%   'XYZ': a 3 column matrix, the coordinate of each gridpt (this
%          could be easily recalculated from IJK, relpt and DX; its
%          mainly kept for debugging or convenience).
%
%   This code relies on a mex-based helper routine that may need to
%   be compiled for your system.

% TODO: write a wrapper or modify this to allow loading directly
% from a ply file.
%   [IJK,DIST,CP,XYZ] = tri2cp(plyname, dx, relpt, p, fd_stenrad)

%addpath('../../surfaces/tri_pig');
%addpath('../../surfaces/readply');

if (nargin < 4)
  relpt = 0;
end
if (nargin < 5)
  disp('tri2cp: defaulting to p=3 and stencil radius 1');
  p = 3;  % interpolation degree
  fd_stenrad = 1;  % Finite difference stencil radius
  calcbw = 1;
elseif (nargin < 6)
  bw = p;
  calcbw = 0;
else
  calcbw = 1;
end

if (calcbw)
  disp('tri2cp: calculuting BW');
  dim = 3;
  % The formula for bw is found in [Ruuth & Merriman 2008] and the
  % 1.00001 is a safety factor.
  bw = 1.00001*sqrt((dim-1)*((p+1)/2)^2 + ((fd_stenrad+(p+1)/2)^2));
  bw = bw * dx;
end

DEBUG_LEVEL = 10;

if (1==0)
  % might be nice to support passing the name of a plyfile instead
  % of V, F
  %PlyFile = 'pig_loop2.ply';
  [Faces, Vertices] = plyread(PlyFile, 'tri');
end

%mex tri2cp_helper.c
%mex -O CFLAGS='\$CFLAGS -Wall' tri2cp_helper.c

trim = 1;

if (all(size(relpt) == [1 1]))
  relpt = [relpt relpt relpt];
end

tic
[IJK,DD,CP,XYZ,CPFACE] = tri2cp_helper(dx, relpt, bw, ...
                                           Faces, Vertices, ...
                                           trim, DEBUG_LEVEL);

% Note: DD is squared distance
DIST = sqrt(DD);

num_grid_points = length(DD)
toc

% a test
assert( max(abs( relpt(1) + (IJK(:,1)-1)*dx - XYZ(:,1) )) == 0)
