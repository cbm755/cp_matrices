function [ijk,dd,cp,xyz,cpface] = tri2cp_helper(dx, relpt, bw, F, V, trim, dbglvl)
%TRI2CP_HELPER  Helper function for TRI2CP
%   Don't call this directly, see TRI2CP instead
%
%   [ijk,dd,cp,xyz,cpface] = helper_tri2cp(dx, relpt, bandwidth, ...
%                            Faces, Vertices, ...
%                            trim, DebugMsgLevel);
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
%   'trim': remove any extra points that light outside 'bandwidth'.
%           Typically, tri2cp will find many more points outside of
%           the specified bandwidth.  If 'trim' is nonzero, these
%           extra points will be discarded.
%
%   'DebugMsgLevel': debug messages from the mex code with a value
%                    less than this will be displayed.
%


%mex tri2cp_helper.c
%mex -O CFLAGS='\$CFLAGS -Wall' tri2cp_helper.c

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
