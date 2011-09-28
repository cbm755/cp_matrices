function [ijk,dd,cp,xyz] = tri2cp_helper(dx, relpt, bw, F, V, gridhashsz, dbglvl)
%TRI2CP_HELPER  Helper function for TRI2CP
%   Don't call this directly, see TRI2CP instead
%
%   [ijk,dd,cp,xyz] = helper_tri2cp(dx, relpt, bandwidth, ...
%                            Faces, Vertices, ...
%                            [expectedgridsz, hashtablesz], ...
%                            DebugMsgLevel);
%
%   'dx': the grid size.
%
%   'relpt': a 3-vector specifying a point from which the integer grid
%            is defined. (if you are debugging and using the full 3D
%            matrix support, then relpt must be lower left corner of
%            data.)
%
%   'bandwidth': the bandwidth that the code should use.
%
%   'Faces', 'Vertices': a triangulation (e.g., see triplot).
%
%   'DebugMsgLevel': debug messages from the mex code with a value
%                    less than this will be displayed.
%
%   TODO: documentation incomplete.



%mex helper_tri2cp.c
%mex -O CFLAGS='\$CFLAGS -Wall' helper_tri2cp.c
%mex -O CFLAGS='\$CFLAGS -Wall -std\=c99' helper_tri2cp.c
disp(' ');
disp('******************************************************');
disp('* TRI2CP_HELPER: compile this code with:');
disp('*');
disp('*      mex tri2cp_helper.c');
disp('*');
disp('* (at the matlab prompt)');
disp('******************************************************');
disp(' ');
error('Compiled code not found');
